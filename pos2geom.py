import linecache

# print lattice constant
for line in range(3,6):
    print ('lattice_vector'+ linecache.getline('./POSCAR',line),end='')

element =linecache.getline('./POSCAR',6).split()
element_numbers =list(map(int,linecache.getline('./POSCAR',7).split()))
coordinate =linecache.getline('./POSCAR',8).strip().upper()

# POSCAR format: Cartesian coordinates
if (coordinate[0] == 'C'):
    for n_element in range(len(element)):
        coor_start = sum(element_numbers[:n_element],9)
        coor_end = sum(element_numbers[:n_element+1],9)
        for geom_line in range(coor_start,coor_end):
            print ('atom' + linecache.getline('./POSCAR',geom_line).rstrip('\n'),element[n_element])

# POSCAR format: Cartesian coordinates
elif (coordinate[0] == 'D'):
    for n_element in range(len(element)):
        coor_start = sum(element_numbers[:n_element],9)
        coor_end = sum(element_numbers[:n_element+1],9)
        for geom_line in range(coor_start,coor_end):
            a = list(map(float,linecache.getline('./POSCAR',3).split()))        
            b = list(map(float,linecache.getline('./POSCAR',4).split()))        
            c = list(map(float,linecache.getline('./POSCAR',5).split()))        

            coor = list(map(float,linecache.getline('./POSCAR',geom_line).split()))

            x = a[0]*coor[0]+b[0]*coor[1]+c[0]*coor[2]
            y = a[1]*coor[0]+b[1]*coor[1]+c[1]*coor[2]
            z = a[2]*coor[0]+b[2]*coor[1]+c[2]*coor[2]
            print ('atom' ,"%15.8f %15.8f %15.8f"%(x,y,z),' ',element[n_element])

#if (coordinate[0] != 'C' and coordinate[0] != 'D'):
else:
    print("The POSCAR is ERROR!!!")
