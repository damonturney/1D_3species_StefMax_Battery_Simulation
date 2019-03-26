"""
Started on Mon Dec 31 2018

Porous electrode theory for Al:EMImCl battery with carbon powder cathode intercalating AlCl4-

This file executes simulations of the solid phase intercalation of a typical pcl
at each locaiton in the electrode, and also a simulation of the diffusion of ions
in the electrolyte.



Simulation of the solid-phase intercalation eqns --- see Damon's notes in Zhang2000
the vector's row/matrix format is setup as shwon below.
the direction of increasing  is into the particles
the column direction is along the electrode
          x_pcl = 0   pppppppppppppppppppppppppppppmmmmmmmmmmmmmmmmmmm x_pcl = 0
                      cccccccccccccccccccccccccccccbbbbbbbbbbbbbbbbbbb
                      lllllllllllllllllllllllllllllrrrrrrrrrrrrrrrrrrr
           x_pcl = 1  12345678910111213141516171819nnnnnnnnnnnnnnnnnnn x_pcl = 1
  current collector | -------------------electrolyte------------------|
        x_elctrlt = 0                                             x_elctrlt = 1


    location 0.0   0.05   0.1  0.15  0.2   0.25  0.3  0.35  0.4   0.45  0.5   0.55  0.6   0.65   0.7  0.75   0.8  0.85   0.9  0.95   1.0  mm
  nodes index 1     2     3     4     5     6     7     8     9     10   11    12    13    14    15    16    17    18    19    20    21  nodes
      boundary|-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----|boundary
     material |carbon pcls carbon pcls carbon pcls carbon pcls carbon pcls carbon pcls carbon pcls carbon  | membrane membrane membran| material

node 1 is at the current collector interface and is infinitesimally inside battery carbon paste (on the lhs) or separator (on the rhs)

For 1-D linear diffusion equation (dm/dt = -D d2m/dx_pcl2) the stability criteria is dt < dx_pcl^2/D

@author: damon
"""


cd("/Users/damon/Desktop/BACKED_UP/WorkFiles/ReportsPublicationsPatents/20181015_JeffXu_ALEMIm_CyclingData_paper")
using SpecialFunctions
using Interpolations
using Plots


#Molar Volumes   A temporary approximation will enforce no net volume change: 4 mv_Al2Cl7 = 7 mv_AlCl4
mv_AlCl4 = 57.14286E-6;   #partial molar volume of  AlCl4- (m3 / mol)
mv_Al2Cl7= 100E-6;   #partial molar volume of Al2Cl7= (m3 / mol)
mv_EMIm  = 150E-6;  #partial molar volume of  EMIm+ (m3 / mol)

#Stefan Maxwell Diffusivities
D10=4.E-7 #in units of m^2/s
D12=4.E-7 #in units of m^2/s
D20=4.E-7 #in units of m^2/s

#Initial Concentrations in the Electrolyte Liquid
elctrlt_thickness=  0.001; # dimensional thicness of the electrolyte domain (meters).   Aluminum metal holds 8 mAh/cm2 in a 10 micron thick film
elctrlt_num_nodes=  31; #number of nodes in the electrolyte domain.         Aluminum metal holds 8 mAh/cm2 in a 10 micron thick film
dx_elctrlt= elctrlt_thickness/(elctrlt_num_nodes-1);  # minus one because one node is on the boundary
separator_thickness=20; #in units of nodes
x_nodes_elctrlt = Array(0:elctrlt_num_nodes-1)*dx_elctrlt;  # locations of nodes 0 to 0.1 cm
num_carbon_nodes=length(x_nodes_elctrlt)-separator_thickness-1
c_liq_AlCl4=    ones(length(x_nodes_elctrlt))*3.0*0.9*1000; #(mol / m3) AlCl4- electrolyte .   Stanford Lin2015 says AlCl4 to EMIm mole ratio varied from 1.1 to 1.8.  The 0.9 is included to make total_molecular_volume_per_volume_of_space be 1 L / 1 L.
c_liq_Al2Cl7=   ones(length(x_nodes_elctrlt))*2.0*0.9*1000; #(mol / m3) Al2Cl7= electrolyte .  Stanford Lin2015 says AlCl4 to EMIm mole ratio varied from 1.1 to 1.8   The 0.9 is included to make total_molecular_volume_per_volume_of_space be 1 L / 1 L.
c_liq_EMIm =    ones(length(x_nodes_elctrlt))*5.0*0.9*1000; #(mol / m3) EMIm+ electrolyte .    Stanford Lin2015 says AlCl4 to EMIm mole ratio varied from 1.1 to 1.8   The 0.9 is included to make total_molecular_volume_per_volume_of_space be 1 L / 1 L.
c_AlCl4_intrcltd=zeros(length(x_nodes_elctrlt));
c_Al_pltd  =    zeros(length(x_nodes_elctrlt));
total_moles_per_liq_volume=c_liq_AlCl4+c_liq_Al2Cl7+c_liq_EMIm;   #(mol/m3)
x_AlCl4=c_liq_AlCl4./total_moles_per_liq_volume;
x_Al2Cl7=c_liq_Al2Cl7./total_moles_per_liq_volume;
x_EMIm=c_liq_EMIm./total_moles_per_liq_volume;
average_molar_volume = x_AlCl4*mv_AlCl4+x_Al2Cl7*mv_Al2Cl7+x_EMIm*mv_EMIm;# m^3 per mole
reaction_stoich_AlCl4 = -1.0*ones(length(x_nodes_elctrlt));
reaction_stoich_AlCl4[end-separator_thickness:end] .= 0.0;
reaction_stoich_AlCl4[end] =  -7.0/3;    # 4Al2Cl7 +3e- -> 7AlCl4 + Al (plated)     divide by 3 because Al3+ requires 3 electrons to plate
reaction_stoich_Al2Cl7= 0.0*ones(length(x_nodes_elctrlt));
reaction_stoich_Al2Cl7[end] = +4.0/3;    # 4Al2Cl7 +3e- -> 7AlCl4 + Al (plated)     divide by 3 because Al3+ requires 3 electrons to plate
reaction_stoich_EMIm  = zeros(length(x_nodes_elctrlt));
reaction_stoich_AlCl4_intrcltd=1.0*ones(length(x_nodes_elctrlt));
reaction_stoich_AlCl4_intrcltd[end-separator_thickness:end] .= 0.0;
reaction_stoich_Al_pltd=zeros(length(x_nodes_elctrlt));
reaction_stoich_Al_pltd[end]=-1.0/3;  # 4Al2Cl7 +3e- -> 7AlCl4 + Al (plated)     divide by 3 because Al3+ requires 3 electrons to plate
z_EMIm=     +1.0;     #charge of the EMIm    (dimensionless)
z_AlCl4=    -1.0;     #charge of the intercalating ion (AlCl4-)
z_Al2Cl7=   -1.0;     #charge of the other ion (Al2Cl7-)

#Designate initial velocities and arrays describing velocity
velocity_AlCl4 = x_nodes_elctrlt*0;  #m/s
velocity_Al2Cl7= x_nodes_elctrlt*0;  #m/s
velocity_EMIm  = x_nodes_elctrlt*0;  #m/s
v_bulk         = x_nodes_elctrlt*0;  #m/s
d_i_sol_dx=      x_nodes_elctrlt*0; #i_sol and i_boundaries have units of Amps per m^2 ;
i_sol =          x_nodes_elctrlt*0;
i_sol_check =    x_nodes_elctrlt*0;
Voltage =        x_nodes_elctrlt*0;
dVdx =           x_nodes_elctrlt*0;
Pressure =       x_nodes_elctrlt*0;
dPdx =           x_nodes_elctrlt*0;
dc_liq_EMIm_dx=  x_nodes_elctrlt*0;
dc_liq_AlCl4_dx= x_nodes_elctrlt*0;
dc_liq_Al2Cl7_dx= x_nodes_elctrlt*0;

#Assorted Variables
T =    298.0;    # Kelvin
F =    96500.0;  # Faradays Constant  (Coulombs/mole)
R =    8.3;      #(J/mole/Kelvin)
Beta = 0.5;      #Butler-Volmer transfer coefficient
r0 =   15E-6    #(m)  radius of cubical particle
volume_fraction_carbon=0.5;
carbon_surface_area_density=6/(2*r0)*volume_fraction_carbon*0.5;  #Units are m2/m3.  The final 0.5 is because some particles will touch each other or not participate for some reason
c_pcl_t0 =   2.28*1000;     #(mol / m3) initial concentration of intercalated-Li throughout pcl
c_host_sites =   24*1000;       #(mol / m3) total concentration of intercalation sites (occupied or unoccpied)
Ustart=Ueq(c_pcl_t0[1]/c_host_sites[1],c_liq_AlCl4[1],R,T,F);# the initial applied potential
Utop=  2.2;      # CV limit
Ubottom=1.8;     # CV limit
kb =   0.05/1000;     ### moles/sec/area/molarity mole/m3   (end result is mole/m2/s)  ...  it has an "kinetic magical" unit of 1 mole/m3 tacked onto it -- see your derivation in Zhang2000
D_pcl =1.0E-12;   ###  (m^2 /s)
dtau=  0.000001;   #dtau=dt*D/r0/r0  (dimensionless)
tau_nodes=Array(0:5000000)*dtau;####500000
saved_time_spacing = 100000; ####10000
vreal= 0.00001; ###(V /s) scan rate
v=     vreal*(r0*r0/D_pcl);

#Create arrays to hold concentrations inside pcls
#dx_real = 1E-4      # in units of cm.  dx is 100 nm
dx_pcl=       0.05; #  non-dimensional x for the pcl simulation
pcl_x_nodes=  Array(0:20)*dx_pcl;  # 0 to 0.1 mm
#dt=       0.001;
dtau=     dtau;      #It was defined above dtau=dt*D/r0/r0
tau_nodes=tau_nodes; #It was defined above
print("dt is ", dtau*r0*r0/D_pcl, " seconds\n");
print("D dt/dx_pcl^2 is ",dtau/dx_pcl/dx_pcl, " and must be less than 0.5\n");
y=        ones(length(pcl_x_nodes),num_carbon_nodes)*c_pcl_t0/c_host_sites;  # array of y , y=c/c_host_sites
dydx_pcl=     zeros(length(pcl_x_nodes),num_carbon_nodes);  # array of dydx_pcl , y=c/c_host_sites
d2ydx_pcl2=   zeros(length(pcl_x_nodes),num_carbon_nodes);  # array of dy2dx_pcl2 , y=c/c_host_sites
dydtau=       zeros(length(pcl_x_nodes),num_carbon_nodes); #
Eta_carbon=   zeros(num_carbon_nodes);


#STABILITY dt < dx_pcl^2/D  , in my case here dt = 0.001  and  dx_pcl^2=0.00001 and D=1
j=1
c_pcl_interface=y[end,:]*c_host_sites;
Ucollector=Ustart;   #must multiply by 1.0 so that Ucollector isn't a duplicate pointer to the memory location of Ustart
Eta_carbon[:] = Ucollector.-Ueq.(y[end,:],c_liq_AlCl4[1:num_carbon_nodes],R,T,F); #Eta_carbon is overpotential
dydx_pcl[end,:]=r0*z_AlCl4*kb*(c_host_sites.-c_pcl_interface)./c_liq_AlCl4[1:num_carbon_nodes]./c_host_sites/D_pcl.*(exp.(-Beta*R*T/F*Eta_carbon) - exp.((1-Beta)*R*T/F*Eta_carbon) );   ### dydx_pcl is non-dimensional here -- see your derivation in Zhang2000
d2ydx_pcl2[1,:]=(2*y[2,:] - 2*y[1,:]) / dx_pcl / dx_pcl;
d2ydx_pcl2[end,:]=(dydx_pcl[end,:]-dydx_pcl[end-1,:])/dx_pcl;
dydtau = d2ydx_pcl2; #+ 2/pcl_x_nodes*dydx_pcl
#print(j, Ucollector,Eta_carbon,c_pcl_interface,dydx_pcl[-1],d2ydx_pcl2[-1])
d_i_sol_dx[1:num_carbon_nodes]=-F*z_AlCl4*D_pcl*dydx_pcl[end,:]*c_host_sites/r0;
d_i_sol_dx[num_carbon_nodes+1:end].=0.0;
i_sol[1] = 0;
for i in 2:length(d_i_sol_dx)-1
    i_sol[i] = i_sol[i-1]+(d_i_sol_dx[i-1]/2+d_i_sol_dx[i]/2)*dx_elctrlt;
end
i_sol[end] = 0;
d_i_sol_dx[end]=(i_sol[end]-i_sol[end-1])/dx_elctrlt;

#Create the arrays that will save the chronological results to disk, and viewed later for insights
saved_time_spacing = saved_time_spacing;  #It was defined above
tau_nodes_saved = cat([1.,2.,3.,4.],saved_time_spacing:saved_time_spacing:length(tau_nodes),dims=1);
y_saved=zeros(length(tau_nodes_saved),length(pcl_x_nodes),num_carbon_nodes);
i_sol_saved=zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
i_sol_check_saved=zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
d_i_sol_dx_saved=zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
Ucollector_saved=zeros((length(tau_nodes_saved)));
y_saved[1,:,:]=y;
tau_saved=zeros(length(tau_nodes_saved));
velocity_AlCl4_saved = zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
velocity_Al2Cl7_saved= zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
velocity_EMIm_saved  = zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
v_bulk_check_saved  = zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
c_liq_AlCl4_saved  =   zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
c_liq_Al2Cl7_saved =   zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
c_liq_EMIm_saved   =   zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
c_AlCl4_intrcltd_saved=zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
c_Al_pltd_saved       =zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
total_molecular_volume_per_volume_of_space=zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
total_molecular_volume_per_volume_of_space[1,:]=(c_liq_AlCl4*mv_AlCl4+c_liq_Al2Cl7*mv_Al2Cl7+c_liq_EMIm*mv_EMIm);  #volume of all the molecules in each node, in units of volume of molecules per volume of space
Eta_carbon_saved   =   zeros((length(tau_nodes_saved)),num_carbon_nodes);
Voltage_saved   =      zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
Pressure_saved =       zeros((length(tau_nodes_saved)),length(x_nodes_elctrlt));
k=1; #for the incrementing index of the saved arrays
y_saved[k,:,:]=y;
i_sol_saved[k,:]=i_sol;
d_i_sol_dx_saved[k,:]=d_i_sol_dx;
Ucollector_saved[k]=Ucollector;
velocity_AlCl4_saved[k,:] =velocity_AlCl4;
velocity_Al2Cl7_saved[k,:]=velocity_Al2Cl7;
velocity_EMIm_saved[k,:]  =velocity_EMIm;
c_liq_AlCl4_saved[k,:]  =  c_liq_AlCl4;
c_liq_Al2Cl7_saved[k,:] =  c_liq_Al2Cl7;
c_liq_EMIm_saved[k,:]   =  c_liq_EMIm;
i_sol_saved[1,:]=i_sol;
j_saved=zeros((length(tau_nodes_saved)));

Ucollector_saved[1]=Ucollector;
println(Ucollector)
for j in 2:length(tau_nodes)   ###loop over non-dimensional time dtau=dt*D/r0/r0
    global Ucollector, v, y, dydtau, c_liq_AlCl4,c_liq_Al2Cl7,c_liq_EMIm,x_EMIm,x_Al2Cl7,x_AlCl4,z_EMIm,z_AlCl4,z_Al2Cl7;
    global k, c_AlCl4_intrcltd, c_Al_pltd, total_moles_per_liq_volume, v_bulk, v_bulk_check, i_sol, i_sol_check;
    #println(j)
    Ucollector=Ucollector+v*dtau;
    if Ucollector >= Utop && sign(v)==1
         v=-v;
        Ucollector=Ucollector+v*dtau*2;
    end
    if Ucollector <= Ubottom && sign(v)==-1
         v=-v;
        Ucollector=Ucollector+v*dtau*2;
    end
    y=y.+dydtau*dtau;       ##### THIS dydtau*dtau IS TOO BIG FOR UNKOWN REASONS!!  COMPARE WITH THE 1D CASE!             #update all of the pcls' solid-phase concentrations, by using the dy/dtau calculated in the previous timestep
    c_pcl_interface=y[end,:]*c_host_sites;      #calculate the surface concentraction of intercalated ions
    Eta_carbon[:] = Ucollector.-Ueq.(y[end,:],c_liq_AlCl4[1:num_carbon_nodes],R,T,F);
    for i in 2:length(pcl_x_nodes)-1;
        dydx_pcl[i,:] = (y[i+1,:] - y[i-1,:]) / 2 / dx_pcl;
    end
    dydx_pcl[end,:]=r0*z_AlCl4*kb*(c_host_sites.-c_pcl_interface)./c_liq_AlCl4[1:num_carbon_nodes]./c_host_sites/D_pcl.*(exp.(-Beta*R*T/F*Eta_carbon) - exp.((1-Beta)*R*T/F*Eta_carbon) );  ### dydx_pcl is non-dimensional here -- see your derivation in Zhang2000
    dydx_pcl[1,:].=0;
    for i in 2:length(pcl_x_nodes)-1
        d2ydx_pcl2[i,:] = (y[i-1,:] - 2*y[i,:] + y[i+1,:]) / dx_pcl / dx_pcl;
    end
    d2ydx_pcl2[1,:]=(2*y[2,:] - 2*y[1,:]) / dx_pcl / dx_pcl;  ## dx_pcl is non-dimensionalized by r0 (pcl cube dimension)
    d2ydx_pcl2[end,:]=(dydx_pcl[end,:]-dydx_pcl[end-1,:])/dx_pcl;
    dydtau = d2ydx_pcl2; #+ 2/pcl_x_nodes*dydx_pcl

    #Calculate the current in the electrolyte and amount of volume/velocity created by the chemical reaction
    d_i_sol_dx[1:num_carbon_nodes]=-F*z_AlCl4*D_pcl*dydx_pcl[end,:]*c_host_sites/r0*carbon_surface_area_density;  #Units of this eqn for d_i_sol_dx work out to C/m3/sec.  If you multiply it by dx it calculates the additional "source" A/m2 to add the the existing A/m2 .
    d_i_sol_dx[num_carbon_nodes+1:end].=0.0;
    i_sol[1] = 0.0;
    for i in 2:length(d_i_sol_dx)
        i_sol[i] = i_sol[i-1]+(d_i_sol_dx[i-1]/2+d_i_sol_dx[i]/2)*dx_elctrlt;
    end
    i_sol[end] = 0;
    d_i_sol_dx[end]= -sum(d_i_sol_dx);

    #Calculate velocities and E and P with the Maxwell-Stefan equations
    velocity_EMIm[1]=   0;
    velocity_AlCl4[1]=  0;
    velocity_Al2Cl7[1]= 0;
    Voltage[1]=0;
    Pressure[1]=0;
    dc_liq_EMIm_dx[1:end]=gradient(c_liq_EMIm[1:end])./dx_elctrlt;
    dc_liq_AlCl4_dx[1:end]=gradient(c_liq_AlCl4[1:end])./dx_elctrlt;
    dc_liq_Al2Cl7_dx[1:end]=gradient(c_liq_Al2Cl7[1:end])./dx_elctrlt;
    for i in 2:length(x_nodes_elctrlt)
        dVdx[i]=             calculate_dVdx(           c_liq_EMIm[i],c_liq_AlCl4[i],c_liq_Al2Cl7[i],dc_liq_EMIm_dx[i],dc_liq_AlCl4_dx[i],dc_liq_Al2Cl7_dx[i],dx_elctrlt,D10,D12,D20,mv_EMIm,mv_AlCl4,mv_Al2Cl7,v_bulk[i],i_sol[i],F,R,T)
        dPdx[i]=             calculate_dPdx(           c_liq_EMIm[i],c_liq_AlCl4[i],c_liq_Al2Cl7[i],dc_liq_EMIm_dx[i],dc_liq_AlCl4_dx[i],dc_liq_Al2Cl7_dx[i],dx_elctrlt,D10,D12,D20,mv_EMIm,mv_AlCl4,mv_Al2Cl7,v_bulk[i],i_sol[i],F,R,T,dVdx[i])
        velocity_Al2Cl7[i]= calculate_velocity_Al2Cl7( c_liq_EMIm[i],c_liq_AlCl4[i],c_liq_Al2Cl7[i],dc_liq_EMIm_dx[i],dc_liq_AlCl4_dx[i],dc_liq_Al2Cl7_dx[i],dx_elctrlt,D10,D12,D20,mv_EMIm,mv_AlCl4,mv_Al2Cl7,v_bulk[i],i_sol[i],F,R,T,dVdx[i],dPdx[i] );
        velocity_AlCl4[i]=  calculate_velocity_AlCl4(  c_liq_EMIm[i],c_liq_AlCl4[i],c_liq_Al2Cl7[i],dc_liq_EMIm_dx[i],dc_liq_AlCl4_dx[i],dc_liq_Al2Cl7_dx[i],dx_elctrlt,D10,D12,D20,mv_EMIm,mv_AlCl4,mv_Al2Cl7,v_bulk[i],i_sol[i],F,R,T,                  velocity_Al2Cl7[i],                  z_EMIm,z_AlCl4,z_Al2Cl7);
        velocity_EMIm[i]=   calculate_velocity_EMIm(   c_liq_EMIm[i],c_liq_AlCl4[i],c_liq_Al2Cl7[i],dc_liq_EMIm_dx[i],dc_liq_AlCl4_dx[i],dc_liq_Al2Cl7_dx[i],dx_elctrlt,D10,D12,D20,mv_EMIm,mv_AlCl4,mv_Al2Cl7,v_bulk[i],i_sol[i],F,R,T,                  velocity_Al2Cl7[i],velocity_AlCl4[i],z_EMIm,z_AlCl4,z_Al2Cl7 );
    end
    #velocity_EMIm[end]=0;
    #velocity_AlCl4[end]=0;
    #velocity_Al2Cl7[end]=0;
    #dVdx[end]=calculate_dVdx(c_liq_EMIm[end],c_liq_AlCl4[end],c_liq_Al2Cl7[end],dc_liq_EMIm_dx[end],dc_liq_AlCl4_dx[end],dc_liq_Al2Cl7_dx[end],dx_elctrlt,D10,D12,D20,mv_EMIm,mv_AlCl4,mv_Al2Cl7,v_bulk[end],i_sol[end],F,R,T)
    #dPdx[end]=calculate_dPdx(c_liq_EMIm[end],c_liq_AlCl4[end],c_liq_Al2Cl7[end],dc_liq_EMIm_dx[end],dc_liq_AlCl4_dx[end],dc_liq_Al2Cl7_dx[end],dx_elctrlt,D10,D12,D20,mv_EMIm,mv_AlCl4,mv_Al2Cl7,v_bulk[end],i_sol[end],F,R,T,dVdx[end])

    #println(velocity_Al2Cl7[25:39])
    v_bulk_check=c_liq_AlCl4.*velocity_AlCl4.*mv_AlCl4.+c_liq_Al2Cl7.*velocity_Al2Cl7*mv_Al2Cl7.+c_liq_EMIm.*velocity_EMIm*mv_EMIm;
    i_sol_check = F*(z_EMIm.*c_liq_EMIm.*velocity_EMIm + z_AlCl4.*c_liq_AlCl4.*velocity_AlCl4 + z_Al2Cl7.*c_liq_Al2Cl7.*velocity_Al2Cl7)

    for i in 2:length(x_nodes_elctrlt)
        Voltage[i] =  Voltage[i-1] +  (dVdx[i-1]+dVdx[i])/2*dx_elctrlt;
        Pressure[i] = Pressure[i-1] + (dPdx[i-1]+dPdx[i])/2*dx_elctrlt;
    end

    ###Update the concentrations  dtau=dt*D/r0/r0
    c_liq_EMIm_old=1.0*c_liq_EMIm;
    c_liq_AlCl4_old=1.0*c_liq_AlCl4;
    c_liq_Al2Cl7_old=1.0*c_liq_Al2Cl7;
    c_AlCl4_intrcltd_old=1.0*c_AlCl4_intrcltd;
    c_Al_pltd_old=1.0*c_Al_pltd;
    c_liq_EMIm=      1.0*update_conc_upwind_advection(dx_elctrlt, dtau*r0*r0/D_pcl, velocity_EMIm,                  c_liq_EMIm_old,        reaction_stoich_EMIm,            d_i_sol_dx, F)
    c_liq_AlCl4=     1.0*update_conc_upwind_advection(dx_elctrlt, dtau*r0*r0/D_pcl, velocity_AlCl4,                 c_liq_AlCl4_old,       reaction_stoich_AlCl4,           d_i_sol_dx, F)
    c_liq_Al2Cl7=    1.0*update_conc_upwind_advection(dx_elctrlt, dtau*r0*r0/D_pcl, velocity_Al2Cl7,                c_liq_Al2Cl7_old,      reaction_stoich_Al2Cl7,          d_i_sol_dx, F)
    c_AlCl4_intrcltd=1.0*update_conc_upwind_advection(dx_elctrlt, dtau*r0*r0/D_pcl, zeros(length(velocity_AlCl4)),  c_AlCl4_intrcltd_old,  reaction_stoich_AlCl4_intrcltd,  d_i_sol_dx, F)
    c_Al_pltd=       1.0*update_conc_upwind_advection(dx_elctrlt, dtau*r0*r0/D_pcl, zeros(length(velocity_AlCl4)),  c_Al_pltd_old,         reaction_stoich_Al_pltd,         d_i_sol_dx, F)

    #println(velocity_Al2Cl7[25:39])

    x_EMIm=c_liq_EMIm/total_moles_per_liq_volume
    x_AlCl4=c_liq_AlCl4/total_moles_per_liq_volume
    x_Al2Cl7=c_liq_Al2Cl7/total_moles_per_liq_volume

    ##### Save data for post-analysis
    if any(j.==tau_nodes_saved)
        k=k+1;
        #print(j)
        #print(" ")
        #println(k)
        println(Ucollector," ",j)#,Eta_carbon,y[1],y[end-2],y[end])
        #println(j, Ucollector,Eta_carbon,c_pcl_interface,dydx_pcl[end],d2ydx_pcl2[end])
        j_saved[k]=j.*1.0;
        tau_saved[k]=tau_nodes[j].*1.0;
        y_saved[k,:,:]=y.*1.0;;
        i_sol_saved[k,:]=i_sol.*1.0;
        i_sol_check_saved[k,:]=i_sol_check.*1.0;
        d_i_sol_dx_saved[k,:]=d_i_sol_dx.*1.0;
        Ucollector_saved[k]=Ucollector.*1.0;
        velocity_AlCl4_saved[k,:] =velocity_AlCl4.*1.0;
        velocity_Al2Cl7_saved[k,:]=velocity_Al2Cl7.*1.0;
        velocity_EMIm_saved[k,:]  =velocity_EMIm.*1.0;
        v_bulk_check_saved[k,:]=v_bulk_check.*1.0;
        Voltage_saved[k,:]=Voltage.*1.0;
        Pressure_saved[k,:]=Pressure.*1.0;
        c_liq_AlCl4_saved[k,:]  =  c_liq_AlCl4.*1.0;
        c_liq_Al2Cl7_saved[k,:] =  c_liq_Al2Cl7.*1.0;
        c_liq_EMIm_saved[k,:]   =  c_liq_EMIm.*1.0;
        c_AlCl4_intrcltd_saved[k,:] = c_AlCl4_intrcltd.*1.0;
        c_Al_pltd_saved[k,:] = c_Al_pltd.*1.0;
        Eta_carbon_saved[k,:]   =  Eta_carbon.*1.0;
        total_molecular_volume_per_volume_of_space[k,:]=c_liq_AlCl4*mv_AlCl4 + c_liq_Al2Cl7*mv_Al2Cl7 + c_liq_EMIm*mv_EMIm;
        total_moles_per_liq_volume=c_liq_EMIm+c_liq_AlCl4+c_liq_Al2Cl7; #
    end
end
print("elapsed time is ",r0*r0/D_pcl*tau_nodes[end], " seconds\n")
plot(Ucollector_saved,i_sol_saved[:,20])
#println(Ucollector_saved)
#println(i_saved[:,1])
plot(x_nodes_elctrlt,i_sol_saved[2,:],title="i_sol")
