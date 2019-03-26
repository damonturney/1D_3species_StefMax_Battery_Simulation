#Plot examples

i=14;
plot(x_nodes_elctrlt, Voltage_saved[i,:],        markershape=:circle,markersize=2,title="Voltage")
plot(x_nodes_elctrlt, Pressure_saved[i,:],       markershape=:circle,markersize=2,title="Pressure")
plot(x_nodes_elctrlt, velocity_EMIm_saved[i,:],  markershape=:circle,markersize=2,title="velocity_EMIm")
plot(x_nodes_elctrlt, velocity_AlCl4_saved[i,:], markershape=:circle,markersize=2,title="velocity_AlCl4")
plot(x_nodes_elctrlt, velocity_Al2Cl7_saved[i,:],markershape=:circle,markersize=2,title="velocity_Al2Cl7")
plot(x_nodes_elctrlt, c_liq_EMIm_saved[i,:],     markershape=:circle,markersize=2,title="c_liq_EMIm",ylims=(4200,4800))
plot(x_nodes_elctrlt, c_liq_AlCl4_saved[i,:],    markershape=:circle,markersize=2,title="c_liq_AlCl4",ylims=(2600,2800))
plot(x_nodes_elctrlt, c_liq_Al2Cl7_saved[i,:],   markershape=:circle,markersize=2,title="c_liq_Al2Cl7",ylims=(1700,1900))
plot(x_nodes_elctrlt, d_i_sol_dx_saved[i,:],     markershape=:circle,markersize=2,title="d_i_sol_dx_saved")
plot(x_nodes_elctrlt, v_bulk_check_saved[i,:],   markershape=:circle,markersize=2,title="v_bulk_check")
plot(x_nodes_elctrlt, i_sol_saved[i,:],          markershape=:circle,markersize=2,title="i_sol")
plot(x_nodes_elctrlt, i_sol_check_saved[i,:],    markershape=:circle,markersize=2,title="i_sol_check")
plot(Ucollector_saved,i_sol_saved[:,20],markershape=:circle,markersize=2,title="I vs V")
#end


plot(Ucollector_saved,Voltage_saved[:,20],markershape=:circle,markersize=2,title="Voltage")
plot(Ucollector_saved,Pressure_saved[:,20],markershape=:circle,markersize=2,title="Pressure")


#y limits    plot(x,y,ylims=(bottom top))
plot(x_nodes_elctrlt,dc_liq_EMIm_dx,markershape=:circle,markersize=2,title="dc_liq_EMIm_dx")
