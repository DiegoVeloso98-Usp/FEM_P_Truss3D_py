import re
import math
import numpy as np
import matplotlib.pyplot as plt

#### Opening data file and storing all numbers in a list "data" ####

file_name = "Est-VonMises.txt"
data = []
n_input = 9
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
    analysed_node = int(input[7])
    analysed_direction = int(input[8])
    
    
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
    

except FileNotFoundError:
    print("Not Found")
except Exception as e:
    print("An error occurred:", e)


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
        S=E*(E/(2*Green+1))
        return S
    else:
        raise ValueError("Invalid deformation type")


    ### Function property returns the vector of properties of the element ###
def Prop(element):
    type=int(element[2])-1                                  #### Minus 1 to match the index of the properties list ####
    prop =(properties[type])     
    return prop



            ### MASS MATRIX ###

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

# Define the convergence criteria
tolerance = 0.000001
max_iterations = 100

# Initialize the displacement vector

Y0=X  
vector_u=np.zeros(npassos)
vector_f_int=np.zeros(npassos)
n_interações=0

        #### Hessian matrix ####
H0 = np.zeros((n_nodes * n_GL, n_nodes * n_GL)) 

################################ DYNAMIC ANALYSIS ################################
for i in range(npassos):
    ipc=i
                       #### External loads vector ####                       
    F_ext = np.array([load * ((i+1)/npassos) for load in ex_loads])

    
    while True: #for it in range(max_iterations): 

        H = np.zeros((n_nodes * n_GL, n_nodes * n_GL))                     #### Hessian matrix ####

        n_interações+=1
        F_int=[0]*len(ex_loads)            #### Equals to zero every iteration ####

        ##### CALCULATION OF Fi AND HESSIAN FOR EACH ELEMENT #####

        for el in range(n_elements):
            node_1=int(connectivity[el][0])                                         ## First and second node of element ##
            node_2=int(connectivity[el][1])
                                                                                        ### (X_0 WONT CHANGE) ###
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
                            elast=float(Prop(connectivity[el])[0])
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
        Y0 += delta_Y0

        #print("H",H)
        #print("g",g)
        #print("F",F_int)
        #print("Y0",Y0)
        #print("delta Y0",delta_Y0)
        #print("norma delta Y0",np.linalg.norm(delta_Y0))
        #print("norma X",np.linalg.norm(X))

        error=(np.linalg.norm(delta_Y0))/np.linalg.norm(X)
        #print("error",error)
        #print("n_interações",n_interações)
        if error < tolerance:
            break



    u_Node= Y0[n_GL*(analysed_node-1)+(analysed_direction-1)]- X [n_GL*(analysed_node-1)+(analysed_direction-1)]
    F_int_analised_node = F_int[n_GL*(analysed_node-1)+(analysed_direction-1)]
    vector_f_int[ipc]=F_int_analised_node

    vector_u[ipc]=u_Node

plt.plot(vector_u, vector_f_int)
plt.xlabel('Displacement')
plt.ylabel('Internal Force')
plt.title('Displacement vs Internal Force')
plt.show()

with open('c:\\Users\\diegoveloso\\Desktop\\MEF_Pos\\V_2\\output.txt', 'w') as file:
    for i in range(len(vector_u)):
        file.write(f"{vector_u[i]}, {vector_f_int[i]}\n")