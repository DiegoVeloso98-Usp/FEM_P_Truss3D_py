import re
import math
import numpy as np
import matplotlib.pyplot as plt

#### Opening data file and storing all numbers in a list "data" ####

file_name = "Est-TreliçaCoda(Força).txt"
data = []
n_input = 11
n_properties = 4  # (Young Modulus, Ni, Área)
n_type_properties = 5  # Number of different types of properties

try:
    with open(file_name, "r") as file:
        for line in file:
            if not line.strip().startswith("#"):
                matches = re.findall(r'-?\d+\.?\d*', line)
                data.extend(matches)



    # Sublist from data, containing the initial information of the .txt
    input = data[0:n_input]
    ##
    n_nodes = int(input[0])
    npassos = int(input[1])
    n_elements = int(input[2])
    n_GL = int(input[3])
    n_Restrictions = int(input[4])
    n_F = int(input[5])
    deformation_Type = int(input[6])
    analysed_node_force = int(input[7])
    analysed_direction_force= int(input[8])
    analysed_node_displacement = int(input[9])
    analysed_direction_displacement = int(input[10])
    

    
    # Connectivity Matrix of size nodes x 3
    connectivity = [
        data[n_input + i * 3:n_input + (3 * (i + 1))]
        for i in range(n_elements)
    ]
    #print(connectivity)

    # Nodes position Matrix
    np_start = n_input + 3 * n_elements

    nodes_position = [
        data[np_start + i * n_GL: np_start + (i + 1) * n_GL]
        for i in range(n_nodes)
    ]
    #print(nodes_position)

    np_end = np_start + n_GL * n_nodes



    # Properties
    properties = [
        data[np_end + i * n_properties: np_end + (i + 1) * n_properties]
        for i in range(n_type_properties)
    ]

    npr_end = np_end + n_properties * n_type_properties

    # Nodes restrictions matrix
    nr_start = npr_end

    nodes_restriction = [
        data[nr_start + i * 3:nr_start + 3 * (1 + i)]   ### 3 due to { node ; direction ; value } ###
        for i in range(n_Restrictions)
    ]
    #print(nodes_restriction)

    nr_end = nr_start + 3 * n_Restrictions  # 3 due to the .txt format { node ; direction ; value }

    # Nodes loads Matrix
    nl_start = nr_end
    loads = [
        data[nl_start + 3 * i:nl_start + 3 * (1 + i)]
        for i in range(n_F)
    ]
    nl_end = nl_start + 3 * n_F 
    print('loads',loads)

    presc_disp=[
        data[nl_end + 3 * i:nl_end + 3 * (1 + i)] for i in range(n_F)
    ]
    print('presc_disp',presc_disp)

except FileNotFoundError:
    print("Not Found")
except Exception as e:
    print("An error occurred:", e)
    print(presc_disp)

#### Global external loads vector ####
ex_loads=[]
actual_node=0


for i in range(n_GL * n_nodes):
    actual_node = math.floor(i / n_GL) + 1
    direction = (i % n_GL) + 1
    found = False                                                                   ## Flag to check if a match is found
    for load in loads:
        if actual_node == int(load[0]) and direction == int(load[1]):
            ex_loads.append(float(load[2]))
            found = True
            break
    if not found:
        ex_loads.append(0)

#### Global nodes position vector ####

X=[]
for i in range(n_nodes): 
    a=i
    for i in range(n_GL):
        X.append(float(nodes_position[a][i]))
    

#X = np.array(nodes_position, dtype=float).flatten()
    

#### Function for Green Strain ####
def length(x1,x2):
    if n_GL==3:
        l=math.sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2+(x1[2]-x2[2])**2)
    elif n_GL==2:
        l=math.sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2)
    else:
        raise ValueError("Invalid number of global coordinates")
    return l


def Egreen_strain(x1,x2,y1,y2):                                                     ## x1,x2,y2,y1 : nodes initial and current positions 
    Egreen= (1/2)*(((length(y1,y2))**2)/((length(x1,x2))**2)-1)
    return Egreen

        #### Functions for calculating the stress for each constitutive model ####


def Stress(E,Green):                        ## E is the vector of properties of the element ##
    
    if deformation_Type==0:                 ## S_Saint-Venant Kirchnoff ##           
        E=float(E[0])
        S=(E*Green)
        return S
    
    elif deformation_Type==1:               ## S_Henky ##
        E=float(E[0])
        S=E*(math.sqrt(1+2*Green)-1)
        return S
    
    elif deformation_Type==2:               ## S_Almaci ##
        E=float(E[0])
        S=E*(Green/(2*Green+1))
        return S
    else:
        raise ValueError("Invalid deformation type")
    
    ### Function that returns the tangent elasticity scalar ##

def Elasticity(E,Green):                        ## E is the vector of properties of the element ##

    if deformation_Type==0:                 ## S_Saint-Venant Kirchnoff ##           
        E=float(E[0])
        TG_E=E
        return TG_E
    
    elif deformation_Type==1:               ## S_Henky ##
        E=float(E[0])
        TG_E=E/(math.sqrt(1+2*Green))
        return TG_E
    
    elif deformation_Type==2:               ## S_Almaci ##
        E=float(E[0])
        TG_E=E/((2*Green+1)**2)
        return TG_E
    else:
        raise ValueError("Invalid deformation type")


    ### Function property returns the vector of properties of the element ###
def Prop(element):
    type=int(element[2])-1                                  #### Minus 1 to match the index of the properties list ####
    prop =(properties[type])     
    return prop



            ### MASS MATRIX ###
''''
m_aux=np.zeros((n_nodes * n_GL, n_nodes * n_GL))
for el in connectivity:
    node_1=int(el[0])
    node_2=int(el[1])
    ro=float(Prop(el)[3])
    A=float(Prop(el)[2])
    l_0=float(length(X[n_GL*(int(el[0])-1):n_GL*int(el[0])],X[n_GL*(int(el[1])-1):n_GL*int(el[1])]))
    print("ro, A, l_0  ",ro*A*l_0/2)
    for i in range(n_GL):
        m_aux[n_GL*(node_1-1)+i][n_GL*(node_1-1)+i]+=ro*A*l_0/2
        m_aux[n_GL*(node_2-1)+i][n_GL*(node_2-1)+i]+=ro*A*l_0/2
'''
###########################################################################################################
#################################### 3D ANALYSIS OUTPUT FILE ##############################################
###########################################################################################################
if n_GL==3:
    output = []
    output.append('MEF-P\n')
    output.append('n_nodes n_elements n_lists\n')
    output.append('#\n')
    output.append(f'{n_nodes} {n_elements} {int(3*npassos)}\n')  
    output.append('cx cy cz dx dy dz\n')
    output.append('#\n')
    for i in range(n_nodes):
        output.append(f'{X[n_GL*i]} {X[n_GL*i+1]} {X[n_GL*i+2]} 0.0 0.0 0.0\n')  ## 3D coordinates ##
    output.append('tipo gaprox no_1.... no_npe group\n')
    output.append('#\n')
    for i in range(n_elements):
        output.append(f'1 1 {connectivity[i][0]} {connectivity[i][1]} 0\n')  ## Element connectivity ##

    output.append('lista de resultados\n')
    output.append('nome da lista\n')
    output.append('dx dy dz valor_cores\n')

#############################################################################################
#############################################################################################
#############################################################################################

# Define the convergence criteria
tolerance = 0.000001
max_iterations = 100

# Initialize the displacement vector

Y0=X  
vector_u=np.zeros(npassos)
vector_f_int=np.zeros(npassos)
vector_f_ext=np.zeros(npassos)
n_interações=0



        #### Hessian matrix ####
H0 = np.zeros((n_nodes * n_GL, n_nodes * n_GL)) 

################################ STATIC ANALYSIS ################################
for type in range(3):
    Y0=X.copy()  
    
    deformation_Type=type
    vector_u=np.zeros(npassos)
    vector_f_int=np.zeros(npassos)

    for i in range(npassos):#
        ipc=i
                        #### External loads vector ####                       
        F_ext = np.array([load * ((i+1)/npassos) for load in ex_loads])

        for displacement in presc_disp:
            node = int(displacement[0])
            dir = int(displacement[1])
            val = float(displacement[2])
            if val != 0:
                Y0[n_GL*(node-1)+(dir-1)]+=val/npassos
        print('Y0',Y0)

        #print("F_ext",F_ext)
        
        while True: #for it in range(max_iterations): 

            H = np.zeros((n_nodes * n_GL, n_nodes * n_GL))                     #### Hessian matrix ####

            n_interações+=1
            F_int=[0]*len(ex_loads)            #### Equals to zero every iteration ####

            ##### CALCULATION OF Fi AND HESSIAN FOR EACH ELEMENT #####

            for el in range(n_elements):
                node_1=int(connectivity[el][0])                                         ## First and second node of element ##
                node_2=int(connectivity[el][1])
                #print('X,,Y',X,Y0)                                                                            ### (X_0 WONT CHANGE) ###
                X0_1 = X [n_GL*(node_1-1) : n_GL*node_1]             ## Initial position of node 1 ##
                X0_2 = X [n_GL*(node_2-1) : n_GL*node_2]             ## Initial position of node 2 ##

                                                                                            ### (Y_0 WILL CHANGE) ###
                Y0_1 = Y0[n_GL*(node_1-1) : n_GL*node_1]                                ## Current position of node 1 ##
                Y0_2 = Y0[n_GL*(node_2-1) : n_GL*node_2]                                ## Current position of node 2 ##

                Egreen=Egreen_strain(X0_1,X0_2,Y0_1,Y0_2)

                #######     INTERNAL FORCE CALCULATION     ########

                for ib in range(2):
                    for ia in range(n_GL):
                        dedyba=(-1)**(ib+1)/length(X0_1,X0_2)*(Y0_2[ia]-Y0_1[ia])
                        #print(dedyba)
                        element=connectivity[el]                                        ## Line 'el' of the connectivity matrix 

                        F_int_elem=float(dedyba*Stress(Prop(element),Egreen)*float(Prop(element)[2]))      ## dedyba * S * A  ##Eq. 2.55 ##

                        index=n_GL*(int(connectivity[el][ib])-1)+ia                     ## Index of F_int_global to add F_int_elem ##
                        F_int[index]+=F_int_elem                                 

                #######     HESSIAN MATRIX CALCULATION     #######
                            
                for ib in range(2):
                    for ia in range(n_GL):
                        dedyba=(-1)**(ib+1)/((length(X0_1,X0_2))**2)*(Y0_2[ia]-Y0_1[ia])

                        for ix in range(2):
                            for ig in range(n_GL): 
                                elast=float(Elasticity((Prop(connectivity[el])),Egreen))
                                area=float(Prop(connectivity[el])[2])
                                dedyxg=(-1)**(ix+1)/(((length(X0_1,X0_2))**2))*(Y0_2[ig]-Y0_1[ig])

                                rh_dy=dedyba*dedyxg*(elast)*area*(length(X0_1,X0_2))        ## AQUI É A ELAST TANGENTE ##

                                d2edybaxg=0
                                if ia==ig:  
                                    d2edybaxg=((-1)**(ib+1))*((-1)**(ix+1))/(length(X0_1,X0_2))**2
                                rh_dy_2=d2edybaxg*Stress(Prop(connectivity[el]),Egreen)*area*(length(X0_1,X0_2))  ## AQUI É A ELAST TANGENTE ##

                                line=n_GL*(int(connectivity[el][ib])-1)+ia
                                column=n_GL*(int(connectivity[el][ix])-1)+ig

                                H[line][column] += float(rh_dy + rh_dy_2)
                
            g= ([float(a - b) for a, b in zip(F_int, F_ext)])

            for restriction in nodes_restriction:
                node = int(restriction[0])
                dir = int(restriction[1])
                val = float(restriction[2])

                if val == 1:
                    for i in range(n_nodes * n_GL):                   ## Loop through the size of H ##
                        H[n_GL * (node-1) + (dir-1)][i] = 0.0         ## Loop through the columns, (dir-1) 1,2,3 -> 0,1,2 ##
                        H[i][n_GL*(node-1)+(dir -1)] = 0.0            ## Loop through the lines ##
                    H[n_GL * (node-1) + (dir-1)][ n_GL * (node-1) + (dir-1)] = 1.0
                g[n_GL * (node-1) + (dir-1)] = 0.0

            H = np.array(H)
            delta_Y0 = np.linalg.solve(H, g)
            delta_Y0 *= -1
            #print('delta Y',delta_Y0) 
            Y0 += delta_Y0
            #print('Y0',Y0)  
            

            error=(np.linalg.norm(delta_Y0))/np.linalg.norm(X)

            if error < tolerance:
                break
  

        u_Node= Y0[n_GL*(analysed_node_displacement-1)+(analysed_direction_displacement-1)]- X [n_GL*(analysed_node_displacement-1)+(analysed_direction_displacement-1)]
        #print('X,,Y',X,Y0)
        F_int_analised_node = F_int[n_GL*(analysed_node_force-1)+(analysed_direction_force-1)]




        vector_f_int[ipc]=F_int_analised_node
        vector_u[ipc]=u_Node
        #print("vector u",vector_u)
        if type==0:
            vector_f_int_svk=vector_f_int
            vector_u_svk=vector_u 
        elif type==1:
            vector_f_int_henky=vector_f_int
            vector_u_henky=vector_u
        elif type==2:
            vector_f_int_almaci=vector_f_int
            vector_u_almaci=vector_u
    ######################    OUTPUT FILE     ######################

        output.append('#\n')
        output.append('desl.x\n')
        for i in range(n_nodes):
                            ##      desl.x               desl.y                      desl.z                         desl.x ##
            output.append(f'{Y0[n_GL*i]-X[n_GL*i]} {Y0[n_GL*i+1]-X[n_GL*i+1]} {Y0[n_GL*i+2]-X[n_GL*i+2]} {Y0[n_GL*i]-X[n_GL*i]} \n')
        output.append('#\n')
        output.append('desl.y\n')
        for i in range(n_nodes):
                            ##      desl.x               desl.y                      desl.z                         desl.y ##
            output.append(f'{Y0[n_GL*i]-X[n_GL*i]} {Y0[n_GL*i+1]-X[n_GL*i+1]} {Y0[n_GL*i+2]-X[n_GL*i+2]} {Y0[n_GL*i+1]-X[n_GL*i+1]} \n')
        output.append('#\n')
        output.append('desl.z\n')
        for i in range(n_nodes):
                            ##      desl.x               desl.y                      desl.z                         desl.z ##
            output.append(f'{Y0[n_GL*i]-X[n_GL*i]} {Y0[n_GL*i+1]-X[n_GL*i+1]} {Y0[n_GL*i+2]-X[n_GL*i+2]} {Y0[n_GL*i+2]-X[n_GL*i+2]} \n')
    
    #output.append('#\n')
    #output.append('sigma.x\n')
'''
    for i in range(n_nodes):
        for el in range(n_elements):
            node_1=int(connectivity[el][0])                                         ## First and second node of element ##
            node_2=int(connectivity[el][1])                                                                           
            X0_1 = X [n_GL*(node_1-1) : n_GL*node_1]             
            X0_2 = X [n_GL*(node_2-1) : n_GL*node_2]                                                                                 
            Y0_1 = Y0[n_GL*(node_1-1) : n_GL*node_1]                                
            Y0_2 = Y0[n_GL*(node_2-1) : n_GL*node_2]                                
            Egreen=Egreen_strain(X0_1,X0_2,Y0_1,Y0_2)
            element=connectivity[el]
            s=Stress(Prop(element),Egreen)
                            ##      desl.x               desl.y                      desl.z          sigma.x ##
            output.append(f'{Y0[n_GL*i]-X[n_GL*i]} {Y0[n_GL*i+1]-X[n_GL*i+1]} {Y0[n_GL*i+2]-X[n_GL*i+2]}  \n')
    output.append('#\n')
'''

arquive=open('acad_out.ogl','w')
arquive.writelines(output)
arquive.close()
    ########### END OUTPUT FILE ###########


plt.plot(-vector_u_svk, -vector_f_int_svk, label='S_Saint-Venant Kirchnoff')
plt.plot(-vector_u_henky, -vector_f_int_henky, label='S_Henky')
plt.plot(-vector_u_almaci, -vector_f_int_almaci, label='S_Almaci')
plt.xlabel('Displacement')
plt.ylabel('Internal Force')
plt.title('Displacement vs Internal Force')
plt.legend()
plt.show()

with open('c:\\Users\\diegoveloso\\Desktop\\MEF_Pos\\V_2\\output.txt', 'w') as file:
    for i in range(len(vector_u)):
        file.write(f"{vector_u[i]}, {vector_f_int[i]}\n")

