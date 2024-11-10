# Heuristic search to the capacitated clustering problem (CCP)
This repository includes the source code of the proposed MA algorithm published in an EJOR paper titled with "Heuristic search to the capacitated clustering problem".

The three sets of 133 benchmark instances are from CCPLIB available from http://www.optsicom.es/ccp/. To facilitate the further research, we upload the instances here.

We made comparisons between MA and some state-of-the-art methods from the following related CCP works:

[1] Martínez-Gavara, A., Campos, V., Gallego, M., Laguna, M., & Martí, R. (2015). Tabu search and GRASP for the capacitated clustering problem. Computational Optimization and Applications, 62, 589-607.

[2] Martinez-Gavara, A., Landa-Silva, D., Campos, V., & Marti, R. (2017). Randomized heuristics for the capacitated clustering problem. Information Sciences, 417, 154-168.

[3] Lai, X., & Hao, J. K. (2016). Iterated variable neighborhood search for the capacitated clustering problem. Engineering Applications of Artificial Intelligence, 56, 102-120.

[4] Brimberg, J., Mladenović, N., Todosijević, R., & Urošević, D. (2019). Solving the capacitated clustering problem with variable neighborhood search. Annals of Operations Research, 272, 289-321.

Please cite our work as:

Zhou, Q., Benlic, U., Wu, Q., & Hao, J. K. (2019). Heuristic search to the capacitated clustering problem. European Journal of Operational Research, 273(2), 464-487.

** Instructions to use the source code of MA

*** To compile:

q.zhou$ make

q.zhou$

*** To run:

q.zhou$ ./MA_CCP ./input_file ./output_sol_file ./out_stat_file

(where "where input_file is the instance name, output_sol_file is a file used to store the solution information, output_stat_file stores the running information")

q.zhou$

*** To clean

q.zhou$ make clean

q.zhou$

