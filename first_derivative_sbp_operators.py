def dx(ux, u, nx, dx, order):
    # summation-by-parts finite difference operators for first derivatives du/dx

    
    m = nx-1
    
    # second order accurate case
    if order==2:
        # calculate partial derivatives on the boundaries:[0, m] with one-sided difference operators
        ux[0, :] = (u[1, :] -  u[0, :])/dx
        ux[m, :] = (u[m, :] -  u[m-1, :])/dx
        
        #calculate partial derivatives in the interior:(1:nx-1)
        for j in range(1, m):
            ux[j, :] = (u[j+1, :] -  u[j-1, :])/(2.0*dx)

        
                   
    # fourth order accurate case        
    if order==4:
        ################################################# 
        # calculate partial derivatives on the boundaries:(0,1,2,3, : m-3, m-2, m-1, m)
        # with one-sided difference operators
        
        ux[0,:] = -24./17*u[0,:] + 59./34*u[1, :]  - 4./17*u[2, :] - 3./34*u[3,:]
        ux[1,:] = -1./2*u[0,:] + 1./2*u[2, :] ;
        ux[2,:] = 4./43*u[0,:] - 59./86*u[1, :]  + 59./86*u[3, :] - 4./43*u[4,:]
        ux[3,:] = 3./98*u[0,:] - 59./98*u[2, :]  + 32./49*u[4, :] - 4./49*u[5,:]


        ux[m,:] = 24./17*u[m,:] - 59./34*u[m-1, :]  + 4./17*u[m-2, :] + 3./34*u[m-3,:]
        ux[m-1,:] = 1./2*u[m,:] - 1./2*u[m-2, :] ;
        ux[m-2,:] = -4./43*u[m,:] + 59./86*u[m-1, :]- 59./86*u[m-3, :]+ 4./43*u[m-4,:]
        ux[m-3,:] = -3./98*u[m,:] + 59./98*u[m-2, :]- 32./49*u[m-4, :]+ 4./49*u[m-5,:]
    
                
        #------------------------------------------------------------------------------------------------------------------------------
        
        for i in range(4, m - 3):
            ux[i,:] = 0.083333333333333*u[i-2,:] - 0.666666666666667*u[i-1,:] + 0.666666666666667*u[i+1,:] - 0.083333333333333*u[i+2,:]

        ux[:,:] = ux/dx

    # sixth order accurate case        
    ################################################# 
       
    if order==6:
        # calculate partial derivatives on the boundaries:(0,1,2,3,4,5,6 : m-6, m-5, m-4, m-3, m-2, m-1, m)
        # with one-sided difference operators 
        ux[0,:] = -1.694834962162858*u[0,:] + 2.245634824947698*u[1,:] - 0.055649692295628*u[2,:] - 0.670383570370653*u[3,:]\
            - 0.188774952148393*u[4,:] + 0.552135032829910*u[5,:] - 0.188126680800077*u[6,:]
        
        ux[1,:] = -0.434411786832708*u[0,:] + 0.107043134706685*u[2,:] + 0.420172642668695*u[3,:] + 0.119957288069806*u[4,:]\
            - 0.328691543801578*u[5,:] + 0.122487487014485*u[6,:] - 0.006557221825386*u[7,:]
        
        ux[2,:] = 0.063307644169533*u[0,:] - 0.629491308812471*u[1,:] + 0.809935419586724*u[3,:] - 0.699016381364484*u[4,:]\
            + 0.850345731199969*u[5,:] - 0.509589652965290*u[6,:] + 0.114508548186019*u[7,:]
        
        ux[3,:] = 0.110198643174386*u[0,:] - 0.357041083340051*u[1,:] - 0.117033418681039*u[2,:] + 0.120870009174558*u[4,:]\
            + 0.349168902725368*u[5,:] - 0.104924741749615*u[6,:] - 0.001238311303608*u[7,:]
        
        ux[4,:] = 0.133544619364965*u[0,:] - 0.438678347579289*u[1,:] + 0.434686341173840*u[2,:] - 0.520172867814934*u[3,:]\
            + 0.049912002176267*u[5,:] + 0.504693510958978*u[6,:] - 0.163985258279827*u[7,:]
        
        ux[5,:] = -0.127754693486067*u[0,:] + 0.393149407857401*u[1,:] - 0.172955234680916*u[2,:] - 0.491489487857764*u[3,:]\
            - 0.016325050231672*u[4,:] + 0.428167552785852*u[6,:] - 0.025864364383975*u[7,:] + 0.013071869997141*u[8,:]
        
        ux[6,:] = 0.060008241515128*u[0,:] - 0.201971348965594*u[1,:] + 0.142885356631256*u[2,:] + 0.203603636754774*u[3,:]\
            - 0.227565385120003*u[4,:] - 0.590259111130048*u[5,:] + 0.757462553894374*u[7,:] - 0.162184436527372*u[8,:]\
            + 0.018020492947486*u[9,:]
        
        ux[7,:] = 0.009910488565285*u[1,:] - 0.029429452176588*u[2,:] + 0.002202493355677*u[3,:] + 0.067773581604826*u[4,:]\
            + 0.032681945726690*u[5,:] - 0.694285851935105*u[6,:] + 0.743286642396343*u[8,:] - 0.148657328479269*u[9,:]\
            + 0.016517480942141*u[10,:]

        ux[m-7,:] =-0.016517480942141*u[m-10,:] + 0.148657328479269*u[m-9,:] - 0.743286642396343*u[m-8,:] + 0.694285851935105*u[m-6,:]\
            - 0.032681945726690*u[m-5,:] - 0.067773581604826*u[m-4,:] - 0.002202493355677*u[m-3,:] + 0.029429452176588*u[m-2,:]\
            - 0.009910488565285*u[m-1,:]

        ux[m-6,:] =-0.018020492947486*u[m-9,:] + 0.162184436527372*u[m-8,:] - 0.757462553894374*u[m-7,:] + 0.590259111130048*u[m-5,:]\
            + 0.227565385120003*u[m-4,:] - 0.203603636754774*u[m-3,:] - 0.142885356631256*u[m-2,:] + 0.201971348965594*u[m-1,:]\
            - 0.060008241515128*u[m,:]

        ux[m-5,:] =-0.013071869997141*u[m-8,:] + 0.025864364383975*u[m-7,:] - 0.428167552785852*u[m-6,:] + 0.016325050231672*u[m-4,:]\
            + 0.491489487857764*u[m-3,:] + 0.172955234680916*u[m-2,:] - 0.393149407857401*u[m-1,:] + 0.127754693486067*u[m,:]

        ux[m-4,:] = 0.163985258279827*u[m-7,:] - 0.504693510958978*u[m-6,:] - 0.049912002176267*u[m-5,:] + 0.520172867814934*u[m-3,:]\
            - 0.434686341173840*u[m-2,:] + 0.438678347579289*u[m-1,:] - 0.133544619364965*u[m,:]

        ux[m-3,:] = 0.001238311303608*u[m-7,:] + 0.104924741749615*u[m-6,:] - 0.349168902725368*u[m-5,:] - 0.120870009174558*u[m-4,:]\
            + 0.117033418681039*u[m-2,:] + 0.357041083340051*u[m-1,:] - 0.110198643174386*u[m,:]

        ux[m-2,:] =-0.114508548186019*u[m-7,:] + 0.509589652965290*u[m-6,:] - 0.850345731199969*u[m-5,:] + 0.699016381364484*u[m-4,:]\
            - 0.809935419586724*u[m-3,:] + 0.629491308812471*u[m-1,:] - 0.063307644169533*u[m,:]

        ux[m-1,:] = 0.006557221825386*u[m-7,:] - 0.122487487014485*u[m-6,:] + 0.328691543801578*u[m-5,:] - 0.119957288069806*u[m-4,:]\
            - 0.420172642668695*u[m-3,:] - 0.107043134706685*u[m-2,:] + 0.434411786832708*u[m,:]

        ux[m,:]   = 0.188126680800077*u[m-6,:] - 0.552135032829910*u[m-5,:] + 0.188774952148393*u[m-4,:] + 0.670383570370653*u[m-3,:]\
            + 0.055649692295628*u[m-2,:] - 2.245634824947698*u[m-1,:] + 1.694834962162858*u[m,:]
        
        for i in range(8, m-7):
            ux[i,:] = -0.016666666666667*u[i-3,:] + 0.15*u[i-2,:] - 0.75*u[i-1,:] + 0.75*u[i+1,:] - 0.15*u[i+2,:] + 0.016666666666667*u[i+3,:]
    
        ux[:,:] = ux[:,:]/dx

def dy(uy, u, ny, dy, order):
    # summation-by-parts finite difference operators for first derivatives du/dy
    
    m = ny-1


    if order == 2:
        uy[:,0]= u[:,1]-u[:,0]
        uy[: m] = u[:,m]-u[:,m-1]
        
        for j in range(1, m):
      
            uy[:, j] = 0.5*(u[:, j+1] - u[:, j-1])
   
        uy[:,:] = 1.0/dy*uy[:,:]

    if order == 4:
        
        uy[:, 0] = -24.0/17.0*u[:,0] + 59.0/34.0*u[:,1]  - 4.0/17.0*u[:, 2] - 3.0/34.0*u[:, 3]
        uy[:, 1] = -1.0/2.0*u[:,0] + 1.0/2.0*u[:,2] 
        uy[:, 2] = 4.0/43.0*u[:,0] - 59.0/86.0*u[:,1]  + 59.0/86.0*u[:,3] - 4.0/43.0*u[:,4]
        uy[:, 3] = 3.0/98.0*u[:,0] - 59.0/98.0*u[:,2]  + 32.0/49.0*u[:,4] - 4.0/49.0*u[:,5]

        uy[:, m] = 24.0/17.0*u[:,m] - 59.0/34.0*u[:,m-1]  + 4.0/17.0*u[:,m-2] + 3.0/34.0*u[:,m-3]
        uy[:, m-1] = 1.0/2.0*u[:,m] - 1.0/2.0*u[:,m-2] 
        uy[:, m-2] = -4.0/43.0*u[:,m] + 59.0/86.0*u[:,m-1]  - 59.0/86.0*u[:,m-3] + 4.0/43.0*u[:,m-4]
        uy[:, m-3] = -3.0/98.0*u[:,m] + 59.0/98.0*u[:,m-2]  - 32.0/49.0*u[:,m-4] + 4.0/49.0*u[:,m-5]

        for j in range(4, m-3):

            uy[:, j] = 1.0/12.0*u[:,j-2] - 2.0/3.0*u[:,j-1]  + 2.0/3.0*u[:,j+1] - 1.0/12.0*u[:,j+2]
    

        uy[:,:] = 1.0/dy*uy[:,:]


    if order == 6:
        
        uy[:,0] = -1.694834962162858*u[:,0] + 2.245634824947698*u[:,1] - 0.055649692295628*u[:,2] - 0.670383570370653*u[:,3]\
            - 0.188774952148393*u[:,4] + 0.552135032829910*u[:,5] - 0.188126680800077*u[:,6]
        
        uy[:,1] = -0.434411786832708*u[:,0] + 0.107043134706685*u[:,2] + 0.420172642668695*u[:,3] + 0.119957288069806*u[:,4]\
            - 0.328691543801578*u[:,5] + 0.122487487014485*u[:,6] - 0.006557221825386*u[:,7]
        
        uy[:,2] = 0.063307644169533*u[:,0] - 0.629491308812471*u[:,1] + 0.809935419586724*u[:,3] - 0.699016381364484*u[:,4]\
            + 0.850345731199969*u[:,5] - 0.509589652965290*u[:,6] + 0.114508548186019*u[:,7]
        
        uy[:,3] = 0.110198643174386*u[:,0] - 0.357041083340051*u[:,1] - 0.117033418681039*u[:,2] + 0.120870009174558*u[:,4]\
            + 0.349168902725368*u[:,5] - 0.104924741749615*u[:,6] - 0.001238311303608*u[:,7]
        
        uy[:,4] = 0.133544619364965*u[:,0] - 0.438678347579289*u[:,1] + 0.434686341173840*u[:,2] - 0.520172867814934*u[:,3]\
            + 0.049912002176267*u[:,5] + 0.504693510958978*u[:,6] - 0.163985258279827*u[:,7]
        
        uy[:,5] = -0.127754693486067*u[:,0] + 0.393149407857401*u[:,1] - 0.172955234680916*u[:,2] - 0.491489487857764*u[:,3]\
            - 0.016325050231672*u[:,4] + 0.428167552785852*u[:,6] - 0.025864364383975*u[:,7] + 0.013071869997141*u[:,8]
        
        uy[:,6] = 0.060008241515128*u[:,0] - 0.201971348965594*u[:,1] + 0.142885356631256*u[:,2] + 0.203603636754774*u[:,3]\
            - 0.227565385120003*u[:,4] - 0.590259111130048*u[:,5] + 0.757462553894374*u[:,7] - 0.162184436527372*u[:,8]\
            + 0.018020492947486*u[:,9]
        
        uy[:,7] = 0.009910488565285*u[:,1] - 0.029429452176588*u[:,2] + 0.002202493355677*u[:,3] + 0.067773581604826*u[:,4]\
            + 0.032681945726690*u[:,5] - 0.694285851935105*u[:,6] + 0.743286642396343*u[:,8] - 0.148657328479269*u[:,9]\
            + 0.016517480942141*u[:,10]

        uy[:,m-7] =-0.016517480942141*u[:,m-10] + 0.148657328479269*u[:,m-9] - 0.743286642396343*u[:,m-8] + 0.694285851935105*u[:,m-6]\
            - 0.032681945726690*u[:,m-5] - 0.067773581604826*u[:,m-4] - 0.002202493355677*u[:,m-3] + 0.029429452176588*u[:,m-2]\
            - 0.009910488565285*u[:,m-1]

        uy[:,m-6] =-0.018020492947486*u[:,m-9] + 0.162184436527372*u[:,m-8] - 0.757462553894374*u[:,m-7] + 0.590259111130048*u[:,m-5]\
            + 0.227565385120003*u[:,m-4] - 0.203603636754774*u[:,m-3] - 0.142885356631256*u[:,m-2] + 0.201971348965594*u[:,m-1]\
            - 0.060008241515128*u[:,m]

        uy[:,m-5] =-0.013071869997141*u[:,m-8] + 0.025864364383975*u[:,m-7] - 0.428167552785852*u[:,m-6] + 0.016325050231672*u[:,m-4]\
            + 0.491489487857764*u[:,m-3] + 0.172955234680916*u[:,m-2] - 0.393149407857401*u[:,m-1] + 0.127754693486067*u[:,m]

        uy[:,m-4] = 0.163985258279827*u[:,m-7] - 0.504693510958978*u[:,m-6] - 0.049912002176267*u[:,m-5] + 0.520172867814934*u[:,m-3]\
            - 0.434686341173840*u[:,m-2] + 0.438678347579289*u[:,m-1] - 0.133544619364965*u[:,m]

        uy[:,m-3] = 0.001238311303608*u[:,m-7] + 0.104924741749615*u[:,m-6] - 0.349168902725368*u[:,m-5] - 0.120870009174558*u[:,m-4]\
            + 0.117033418681039*u[:,m-2] + 0.357041083340051*u[:,m-1] - 0.110198643174386*u[:,m]

        uy[:,m-2] =-0.114508548186019*u[:,m-7] + 0.509589652965290*u[:,m-6] - 0.850345731199969*u[:,m-5] + 0.699016381364484*u[:,m-4]\
            - 0.809935419586724*u[:,m-3] + 0.629491308812471*u[:,m-1] - 0.063307644169533*u[:,m]

        uy[:,m-1] = 0.006557221825386*u[:,m-7] - 0.122487487014485*u[:,m-6] + 0.328691543801578*u[:,m-5] - 0.119957288069806*u[:,m-4]\
            - 0.420172642668695*u[:,m-3] - 0.107043134706685*u[:,m-2] + 0.434411786832708*u[:,m]

        uy[:,m]   = 0.188126680800077*u[:,m-6] - 0.552135032829910*u[:,m-5] + 0.188774952148393*u[:,m-4] + 0.670383570370653*u[:,m-3]\
            + 0.055649692295628*u[:,m-2] - 2.245634824947698*u[:,m-1] + 1.694834962162858*u[:,m]
                       
        for j in range(8, m-7):
            uy[:, j] = -0.016666666666666667*u[:,j-3] + 0.15*u[:,j-2] - 0.75*u[:,j-1] + 0.75*u[:,j+1] - 0.15*u[:,j+2] + 0.01666666666666667*u[:,j+3]
            
    
        uy[:,:] = (1.0/dy)*uy[:,:];


        

    
def dx2d(ux, u, nx, i, j, dx, order):
    # summation-by-parts finite difference operators for first derivatives du/dx
    
    m = nx-1
    
    # second order accurate case
    if order==2:
        # calculate partial derivatives on the boundaries:[0, m] with one-sided difference operators

        if i == 0:
            ux[:] = (u[1, j, :] -  u[0,j, :])/dx

        if i == m:
            ux[:] = (u[m,j, :] -  u[m-1,j, :])/dx
        
        #calculate partial derivatives in the interior:(1:nx-1)
        if (i > 0) and (i < m):
            ux[:] = (u[i+1, j,:] -  u[i-1, j,:])/(2.0*dx)

        
                   
    # fourth order accurate case        
    if order==4:
        ################################################# 
        # calculate partial derivatives on the boundaries:(0,1,2,3, : m-3, m-2, m-1, m)
        # with one-sided difference operators

        if i == 0:
            ux[:] = -24./17*u[0,j,:] + 59./34*u[1, j,:]  - 4./17*u[2, j,:] - 3./34*u[3,j,:]

        if i == 1:
            ux[:] = -1./2*u[0,j,:] + 1./2*u[2, j,:]

        if i == 2:
            ux[:] = 4./43*u[0,j,:] - 59./86*u[1, j,:]  + 59./86*u[3, j,:] - 4./43*u[4,j,:]

        if i == 3:    
            ux[:] = 3./98*u[0,j,:] - 59./98*u[2, j,:]  + 32./49*u[4, j,:] - 4./49*u[5,j,:]


        if i == m:
            ux[:] = 24./17*u[m,j,:] - 59./34*u[m-1, j,:]  + 4./17*u[m-2, j,:] + 3./34*u[m-3,j,:]

        if i == m-1:
            ux[:] = 1./2*u[m,j,:] - 1./2*u[m-2, j,:]

        if i == m-2:
            ux[:] = -4./43*u[m,j,:] + 59./86*u[m-1, j,:]- 59./86*u[m-3, j,:]+ 4./43*u[m-4,j,:]

        if i == m-3:
            ux[:] = -3./98*u[m,j,:] + 59./98*u[m-2, j,:]- 32./49*u[m-4, j,:]+ 4./49*u[m-5,j,:]
    
                
        #------------------------------------------------------------------------------------------------------------------------------
        
        if (i > 3) and (i<m-3):
            ux[:] = 0.083333333333333*u[i-2,j,:] - 0.666666666666667*u[i-1,j,:] + 0.666666666666667*u[i+1,j,:] - 0.083333333333333*u[i+2,j,:]

        ux[:] = ux[:]/dx

    # sixth order accurate case        
    ################################################# 
       
    if order==6:
        # calculate partial derivatives on the boundaries:(0,1,2,3,4,5,6 : m-6, m-5, m-4, m-3, m-2, m-1, m)
        # with one-sided difference operators

        if i == 0:
            ux[:] = -1.694834962162858*u[0,j,:] + 2.245634824947698*u[1,j,:] - 0.055649692295628*u[2,j,:] - 0.670383570370653*u[3,j,:]\
                    - 0.188774952148393*u[4,j,:] + 0.552135032829910*u[5,j,:] - 0.188126680800077*u[6,j,:]
            
        if i == 1:
            ux[:] = -0.434411786832708*u[0,j,:] + 0.107043134706685*u[2,j,:] + 0.420172642668695*u[3,j,:] + 0.119957288069806*u[4,j,:]\
                      - 0.328691543801578*u[5,j,:] + 0.122487487014485*u[6,j,:] - 0.006557221825386*u[7,j,:]

        if i == 2:
            ux[:] = 0.063307644169533*u[0,j,:] - 0.629491308812471*u[1,j,:] + 0.809935419586724*u[3,j,:] - 0.699016381364484*u[4,j,:]\
                    + 0.850345731199969*u[5,j,:] - 0.509589652965290*u[6,j,:] + 0.114508548186019*u[7,j,:]

        if i == 3:
            ux[:] = 0.110198643174386*u[0,j,:] - 0.357041083340051*u[1,j,:] - 0.117033418681039*u[2,j,:] + 0.120870009174558*u[4,j,:]\
                    + 0.349168902725368*u[5,j,:] - 0.104924741749615*u[6,j,:] - 0.001238311303608*u[7,j,:]

        if i == 4:
            ux[:] = 0.133544619364965*u[0,j,:] - 0.438678347579289*u[1,j,:] + 0.434686341173840*u[2,j,:] - 0.520172867814934*u[3,j,:]\
                    + 0.049912002176267*u[5,j,:] + 0.504693510958978*u[6,j,:] - 0.163985258279827*u[7,j,:]

        if i == 5:
            ux[:] = -0.127754693486067*u[0,j,:] + 0.393149407857401*u[1,j,:] - 0.172955234680916*u[2,j,:] - 0.491489487857764*u[3,j,:]\
                    - 0.016325050231672*u[4,j,:] + 0.428167552785852*u[6,j,:] - 0.025864364383975*u[7,j,:] + 0.013071869997141*u[8,j,:]

        if i == 6:
            ux[:] = 0.060008241515128*u[0,j,:] - 0.201971348965594*u[1,j,:] + 0.142885356631256*u[2,j,:] + 0.203603636754774*u[3,j,:]\
                    - 0.227565385120003*u[4,j,:] - 0.590259111130048*u[5,j,:] + 0.757462553894374*u[7,j,:] - 0.162184436527372*u[8,j,:]\
                    + 0.018020492947486*u[9,j,:]
        
        if i == 7:
            ux[:] = 0.009910488565285*u[1,j,:] - 0.029429452176588*u[2,j,:] + 0.002202493355677*u[3,j,:] + 0.067773581604826*u[4,j,:]\
                      + 0.032681945726690*u[5,j,:] - 0.694285851935105*u[6,j,:] + 0.743286642396343*u[8,j,:] - 0.148657328479269*u[9,j,:]\
                      + 0.016517480942141*u[10,j,:]

        if i == m-7:
            ux[:] =-0.016517480942141*u[m-10,j,:] + 0.148657328479269*u[m-9,j,:] - 0.743286642396343*u[m-8,j,:] + 0.694285851935105*u[m-6,j,:]\
                - 0.032681945726690*u[m-5,j,:] - 0.067773581604826*u[m-4,j,:] - 0.002202493355677*u[m-3,j,:] + 0.029429452176588*u[m-2,j,:]\
                - 0.009910488565285*u[m-1,j,:]

        if i == m-6:
            ux[:] =-0.018020492947486*u[m-9,j,:] + 0.162184436527372*u[m-8,j,:] - 0.757462553894374*u[m-7,j,:] + 0.590259111130048*u[m-5,j,:]\
                + 0.227565385120003*u[m-4,j,:] - 0.203603636754774*u[m-3,j,:] - 0.142885356631256*u[m-2,j,:] + 0.201971348965594*u[m-1,j,:]\
                - 0.060008241515128*u[m,j,:]

        if i == m-5:
            ux[:] =-0.013071869997141*u[m-8,j,:] + 0.025864364383975*u[m-7,j,:] - 0.428167552785852*u[m-6,j,:] + 0.016325050231672*u[m-4,j,:]\
                + 0.491489487857764*u[m-3,j,:] + 0.172955234680916*u[m-2,j,:] - 0.393149407857401*u[m-1,j,:] + 0.127754693486067*u[m,j,:]

        if i == m-4:
            ux[:] = 0.163985258279827*u[m-7,j,:] - 0.504693510958978*u[m-6,j,:] - 0.049912002176267*u[m-5,j,:] + 0.520172867814934*u[m-3,j,:]\
                    - 0.434686341173840*u[m-2,j,:] + 0.438678347579289*u[m-1,j,:] - 0.133544619364965*u[m,j,:]

        if i == m-3:
            ux[:] = 0.001238311303608*u[m-7,j,:] + 0.104924741749615*u[m-6,j,:] - 0.349168902725368*u[m-5,j,:] - 0.120870009174558*u[m-4,j,:]\
                    + 0.117033418681039*u[m-2,j,:] + 0.357041083340051*u[m-1,j,:] - 0.110198643174386*u[m,j,:]

        if i == m-2:
            ux[:] =-0.114508548186019*u[m-7,j,:] + 0.509589652965290*u[m-6,j,:] - 0.850345731199969*u[m-5,j,:] + 0.699016381364484*u[m-4,j,:]\
                - 0.809935419586724*u[m-3,j,:] + 0.629491308812471*u[m-1,j,:] - 0.063307644169533*u[m,j,:]

        if i == m-1:
            ux[:] = 0.006557221825386*u[m-7,j,:] - 0.122487487014485*u[m-6,j,:] + 0.328691543801578*u[m-5,j,:] - 0.119957288069806*u[m-4,j,:]\
                    - 0.420172642668695*u[m-3,j,:] - 0.107043134706685*u[m-2,j,:] + 0.434411786832708*u[m,j,:]

        if i == m:
            ux[:]   = 0.188126680800077*u[m-6,j,:] - 0.552135032829910*u[m-5,j,:] + 0.188774952148393*u[m-4,j,:] + 0.670383570370653*u[m-3,j,:]\
                      + 0.055649692295628*u[m-2,j,:] - 2.245634824947698*u[m-1,j,:] + 1.694834962162858*u[m,j,:]
        
        if (i > 7) and (i < m-7):
            ux[:] = -0.016666666666667*u[i-3,j,:] + 0.15*u[i-2,j,:] - 0.75*u[i-1,j,:] + 0.75*u[i+1,j,:] - 0.15*u[i+2,j,:] + 0.016666666666667*u[i+3,j,:]
    
        ux[:] = ux[:]/dx

def dy2d(uy, u, ny, i, j, dy, order):
    # summation-by-parts finite difference operators for first derivatives du/dy
    
    m = ny-1


    if order == 2:

        if j == 0:
            uy[:]= u[i,1,:]-u[i,0,:]

        if j == m:    
            uy[:] = u[i,m,:]-u[i,m-1,:]
        
        if (j > 0) and j < m:
      
            uy[:] = 0.5*(u[i, j+1,:] - u[i, j-1,:])
   
        uy[:] = 1.0/dy*uy[:]

    if order == 4:
        
        if j == 0:
            uy[:] = -24.0/17.0*u[i,0,:] + 59.0/34.0*u[i,1,:]  - 4.0/17.0*u[i, 2,:] - 3.0/34.0*u[i, 3,:]
            
        if j == 1:
            uy[:] = -1.0/2.0*u[i,0,:] + 1.0/2.0*u[i,2,:]

        if j == 2:
            uy[:] = 4.0/43.0*u[i,0,:] - 59.0/86.0*u[i,1,:]  + 59.0/86.0*u[i,3,:] - 4.0/43.0*u[i,4,:]

        if j == 3:
             uy[:] = 3.0/98.0*u[i,0,:] - 59.0/98.0*u[i,2,:]  + 32.0/49.0*u[i,4,:] - 4.0/49.0*u[i,5,:]

        if j == m:
            uy[:] = 24.0/17.0*u[i,m,:] - 59.0/34.0*u[i,m-1,:]  + 4.0/17.0*u[i,m-2,:] + 3.0/34.0*u[i,m-3,:]

        if j == m-1:    
            uy[:] = 1.0/2.0*u[i,m,:] - 1.0/2.0*u[i,m-2,:]

        if j == m-2:    
            uy[:] = -4.0/43.0*u[i,m,:] + 59.0/86.0*u[i,m-1,:]  - 59.0/86.0*u[i,m-3,:] + 4.0/43.0*u[i,m-4,:]

        if j == m-3:    
            uy[:] = -3.0/98.0*u[i,m,:] + 59.0/98.0*u[i,m-2,:]  - 32.0/49.0*u[i,m-4,:] + 4.0/49.0*u[i,m-5,:]

        if (j> 3) and (j<m-3):
            uy[:] = 1.0/12.0*u[i,j-2,:] - 2.0/3.0*u[i,j-1,:]  + 2.0/3.0*u[i,j+1,:] - 1.0/12.0*u[i,j+2,:]
    

        uy[:] = 1.0/dy*uy[:]


    if order == 6:

        if j == 0:
            uy[:] = -1.694834962162858*u[i,0,:] + 2.245634824947698*u[i,1,:] - 0.055649692295628*u[i,2,:] - 0.670383570370653*u[i,3,:]\
                    - 0.188774952148393*u[i,4,:] + 0.552135032829910*u[i,5,:] - 0.188126680800077*u[i,6,:]

        if j == 1:
            uy[:] = -0.434411786832708*u[i,0,:] + 0.107043134706685*u[i,2,:] + 0.420172642668695*u[i,3,:] + 0.119957288069806*u[i,4,:]\
                    - 0.328691543801578*u[i,5,:] + 0.122487487014485*u[i,6,:] - 0.006557221825386*u[i,7,:]

        if j == 2:
            uy[:] = 0.063307644169533*u[i,0,:] - 0.629491308812471*u[i,1,:] + 0.809935419586724*u[i,3,:] - 0.699016381364484*u[i,4,:]\
                + 0.850345731199969*u[i,5,:] - 0.509589652965290*u[i,6,:] + 0.114508548186019*u[i,7,:]
        
        if j == 3:
            uy[:] = 0.110198643174386*u[i,0,:] - 0.357041083340051*u[i,1,:] - 0.117033418681039*u[i,2,:] + 0.120870009174558*u[i,4,:]\
                    + 0.349168902725368*u[i,5,:] - 0.104924741749615*u[i,6,:] - 0.001238311303608*u[i,7,:]

        if j == 4:
            uy[:] = 0.133544619364965*u[i,0,:] - 0.438678347579289*u[i,1,:] + 0.434686341173840*u[i,2,:] - 0.520172867814934*u[i,3,:]\
                    + 0.049912002176267*u[i,5,:] + 0.504693510958978*u[i,6,:] - 0.163985258279827*u[i,7,:]

        if j == 5:
            uy[:] = -0.127754693486067*u[i,0,:] + 0.393149407857401*u[i,1,:] - 0.172955234680916*u[i,2,:] - 0.491489487857764*u[i,3,:]\
                      - 0.016325050231672*u[i,4,:] + 0.428167552785852*u[i,6,:] - 0.025864364383975*u[i,7,:] + 0.013071869997141*u[i,8,:]

        if j == 6:
            uy[:] = 0.060008241515128*u[i,0,:] - 0.201971348965594*u[i,1,:] + 0.142885356631256*u[i,2,:] + 0.203603636754774*u[i,3,:]\
                      - 0.227565385120003*u[i,4,:] - 0.590259111130048*u[i,5,:] + 0.757462553894374*u[i,7,:] - 0.162184436527372*u[i,8,:]\
                      + 0.018020492947486*u[i,9,:]

        if j == 7:
            uy[:] = 0.009910488565285*u[i,1,:] - 0.029429452176588*u[i,2,:] + 0.002202493355677*u[i,3,:] + 0.067773581604826*u[i,4,:]\
                      + 0.032681945726690*u[i,5,:] - 0.694285851935105*u[i,6,:] + 0.743286642396343*u[i,8,:] - 0.148657328479269*u[i,9,:]\
                      + 0.016517480942141*u[i,10,:]

        if j == m-7:
            uy[:] =-0.016517480942141*u[i,m-10,:] + 0.148657328479269*u[i,m-9,:] - 0.743286642396343*u[i,m-8,:] + 0.694285851935105*u[i,m-6,:]\
                - 0.032681945726690*u[i,m-5,:] - 0.067773581604826*u[i,m-4,:] - 0.002202493355677*u[i,m-3,:] + 0.029429452176588*u[i,m-2,:]\
                - 0.009910488565285*u[i,m-1,:]

        if j == m-6:
            uy[:] =-0.018020492947486*u[i,m-9,:] + 0.162184436527372*u[i,m-8,:] - 0.757462553894374*u[i,m-7,:] + 0.590259111130048*u[i,m-5,:]\
                + 0.227565385120003*u[i,m-4,:] - 0.203603636754774*u[i,m-3,:] - 0.142885356631256*u[i,m-2,:] + 0.201971348965594*u[i,m-1,:]\
                - 0.060008241515128*u[i,m,:]

        if j == m-5:
            uy[:] =-0.013071869997141*u[i,m-8,:] + 0.025864364383975*u[i,m-7,:] - 0.428167552785852*u[i,m-6,:] + 0.016325050231672*u[i,m-4,:]\
            + 0.491489487857764*u[i,m-3,:] + 0.172955234680916*u[i,m-2,:] - 0.393149407857401*u[i,m-1,:] + 0.127754693486067*u[i,m,:]

        if j == m-4:
            uy[:] = 0.163985258279827*u[i,m-7,:] - 0.504693510958978*u[i,m-6,:] - 0.049912002176267*u[i,m-5,:] + 0.520172867814934*u[i,m-3,:]\
                          - 0.434686341173840*u[i,m-2,:] + 0.438678347579289*u[i,m-1,:] - 0.133544619364965*u[i,m,:]

        if j == m-3:
            uy[:] = 0.001238311303608*u[i,m-7,:] + 0.104924741749615*u[i,m-6,:] - 0.349168902725368*u[i,m-5,:] - 0.120870009174558*u[i,m-4,:]\
                          + 0.117033418681039*u[i,m-2,:] + 0.357041083340051*u[i,m-1,:] - 0.110198643174386*u[i,m,:]

        if j == m-2:
            uy[:] =-0.114508548186019*u[i,m-7,:] + 0.509589652965290*u[i,m-6,:] - 0.850345731199969*u[i,m-5,:] + 0.699016381364484*u[i,m-4,:]\
                - 0.809935419586724*u[i,m-3,:] + 0.629491308812471*u[i,m-1,:] - 0.063307644169533*u[i,m,:]

        if j == m-1:
            uy[:] = 0.006557221825386*u[i,m-7,:] - 0.122487487014485*u[i,m-6,:] + 0.328691543801578*u[i,m-5,:] - 0.119957288069806*u[i,m-4,:]\
                    - 0.420172642668695*u[i,m-3,:] - 0.107043134706685*u[i,m-2,:] + 0.434411786832708*u[i,m,:]
        
        if j == m:
            uy[:]   = 0.188126680800077*u[i,m-6,:] - 0.552135032829910*u[i,m-5,:] + 0.188774952148393*u[i,m-4,:] + 0.670383570370653*u[i,m-3,:]\
                    + 0.055649692295628*u[i,m-2,:] - 2.245634824947698*u[i,m-1,:] + 1.694834962162858*u[i,m,:]
                       
        if (j> 7) and (j<m-7):
            uy[:] = -0.016666666666666667*u[i,j-3,:] + 0.15*u[i,j-2,:] - 0.75*u[i,j-1,:] + 0.75*u[i,j+1,:] - 0.15*u[i,j+2,:] + 0.01666666666666667*u[i,j+3,:]
            
    
        uy[:] = (1.0/dy)*uy[:]

    
