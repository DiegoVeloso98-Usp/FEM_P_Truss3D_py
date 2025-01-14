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
    n_steps = int(input[1])
    n_elements = int(input[2])
    n_GL = int(input[3])
    n_Restrictions = int(input[4])
    n_F = int(input[5])
    deformation_Type = int(input[6])
    analysed_node = int(input[7])
    analysed_direction = int(input[8])
 
    
    # CONNECTIVITY MATRIX (n_nodes) x 3 #
    connectivity = [
        data[n_input + i * 3:n_input + (3 * (i + 1))]
        for i in range(n_elements)
    ]

    # NODES POSITION MATRIX #
    np_start = n_input + 3 * n_elements

    nodes_position = [
        data[np_start + i * n_GL: np_start + (i + 1) * n_GL]
        for i in range(n_nodes)
    ]

    np_end = np_start + n_GL * n_nodes

    # PROPERTIES #
    properties = [
        data[np_end + i * n_properties: np_end + (i + 1) * n_properties]
        for i in range(n_type_properties)
    ]

    npr_end = np_end + n_properties * n_type_properties

    # NODES RESTRICTIONS MATRIX # 
    nr_start = npr_end

    nodes_restriction = [
        data[nr_start + i * 3:nr_start + 3 * (1 + i)]   ### 3 due to { node ; direction ; value } ###
        for i in range(n_Restrictions)
    ]

    nr_end = nr_start + 3 * n_Restrictions  # 3 due to the .txt format { node ; direction ; value }

    # NODES LOADS MATRIX #
    nl_start = nr_end
    loads = [
        data[nl_start + 3 * i:nl_start + 3 * (1 + i)]
        for i in range(n_F)
    ]

except FileNotFoundError:
    exit("File Not Found")
 
except Exception as e:
    exit("An error occurred:", e)
    


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

#X = np.array(nodes_position, dtype=float).flatten()

X=[]
for i in range(n_nodes): 
    a=i
    for i in range(n_GL):
        X.append(float(nodes_position[a][i]))

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


#n=np.ones(n_nodes * n_GL)


# Define the convergence criteria

tolerance = 1e-6
max_iterations = 100

        #### Initializing Vectors ####

vector_u=np.zeros(n_steps)
vector_f_int=np.zeros(n_steps)
n_interações=0


Y=X  
H0 = np.zeros((n_nodes * n_GL, n_nodes * n_GL)) 

F_ext=np.array(ex_loads)

################################ DYNAMIC ANALYSIS ################################

for i in range(n_steps):     ## TIME STEPS ##
    ipc=i
    F_ext = (i+1)/(n_steps)*F_ext                                    #### Eq. 2.100 ####


    while True: # for it in range(max_iterations):
        H = H0                                                                 #### H = 0 every iteration ####
        F_int=[0]*len(ex_loads)                                             #### F_int = 0 every iteration ####
        n_interações+=1

        ##### CALCULATION OF Fi AND HESSIAN FOR EACH ELEMENT #####

        for el in range(n_elements):
            node_1=int(connectivity[el][0])                                         ## First and second node of element ##
            node_2=int(connectivity[el][1])
                                                                                        ### (X_0 WONT CHANGE) ###
            X0_1 = X [n_GL*(node_1-1) : n_GL*node_1]             ## Initial position of node 1 ##
            X0_2 = X [n_GL*(node_2-1) : n_GL*node_2]             ## Initial position of node 2 ##
                                                                                        ### (Y_0 WILL CHANGE) ###
            Y0_1 = Y[n_GL*(node_1-1) : n_GL*node_1]                                ## Current position of node 1 ##
            Y0_2 = Y[n_GL*(node_2-1) : n_GL*node_2]                                ## Current position of node 2 ##

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
                    #print(F_int) #'''PAREI A VERIFICAÇAO AQUI (F_INT ESTÁ DANDO 0 PARA TODOS OS ELEMENTOS)'''

                                                     
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
                            H[line][column] += float(rh_dy + rh_dy_2)   # Static Hessian #

              
        
        g= ([float(a - b) for a, b in zip(F_int, F_ext)])  ## Residual vector ##

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
        
        
        delta_Y = np.linalg.solve(H, g)
        delta_Y *= -1
        Y += delta_Y
        error=(np.linalg.norm(delta_Y))/np.linalg.norm(X)
        if error < tolerance:
            break

    u_Node= Y [n_GL*(analysed_node-1)+(analysed_direction-1)]- X [n_GL*(analysed_node-1)+(analysed_direction-1)]
    F_int_analised_node = F_int[n_GL*(analysed_node-1)+(analysed_direction-1)]
    vector_f_int[ipc]=F_int_analised_node
    vector_u[ipc]=u_Node

plt.plot(vector_u, vector_f_int)
plt.xlabel('Displacement')
plt.ylabel('Internal Force')
plt.title('Displacement vs Internal Force')
plt.show()
with open('G:\Meu Drive\MEF_Pos\V_6 - Static  - Validation\Est-TorreEngastada-output', 'w') as file:
    for i in range(len(vector_u)):
        file.write(f"{vector_u[i]}, {vector_f_int[i]}\n")