import numpy as np
import test_rm_2_n as ld

F = np.array([
    [ [1, 2], [3, 4] ], 
    [ [5, 6], [7, 8] ]  
])
a  = ld.FFT_one_step(F, 1)
print(a[0][0])




a = ld.sums_RM2(np.array([0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1,
0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0,
0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0,
1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0,
1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1,
0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]), 0.3)
result_distance = [i[1] for i in a]
max_value = max(result_distance, key=abs)
print("the maximum value is:"+str(max_value))
index_quadratic_nearest = result_distance.index(max_value)
print("the nearest Boolean function is:")
print(a[index_quadratic_nearest][0].lins,a[index_quadratic_nearest][0].lin)
for i in range(len(a[index_quadratic_nearest][0].lins)):
    if a[index_quadratic_nearest][0].lins[i]:
        print("x"+str(i)+"*("+str(ld.Linear(a[index_quadratic_nearest][0].lins[i]))+")",end="+")
    # print(a[index_quadratic_nearest][0].lins[2], a[index_quadratic_nearest][0].lin)
print(a[index_quadratic_nearest][0].lin)





# tt:= [];
# for x9,x8,x7,x6,x5,x4,x3,x2,x1,x0 in GF(2) do 
#     Append(~tt,x6*(x2+x3+x4+x5)+x7*(x0+x1+x3+x4+x5+x6)+x8*(x0+x1+x3+x7)+x9*(x0+x1+x5+x6+x7)+x1+x3+x4+x6);
# end for;



# tt:=[];
# for x0,x1,x2,x3,x4,x5,x6 in GF(2) do
#     Append(~tt, );
# end for;
