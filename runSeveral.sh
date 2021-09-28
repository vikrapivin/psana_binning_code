#mpiexec -n 6 python process_cube_h5.py --time_low -6 --time_high 51 --laser_off 1 126
#for ii in {136..138}
#do
##echo $ii
##mpiexec -n 6 python process_cube_h5.py --laser_off 1 $ii
#mpiexec -n 6 python spt_process_cube.py --time_low -6 --time_high 51 --laser_off 1 $ii
#done
