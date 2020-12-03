%% Beispiel
G = [.5 , .2 , .5 , .7 , .2 , .5 , .7 , 1;
     .1 , .3 , .3 , .3 , .5 , .5 , .5 , .5];
B = [.28 , .63 , .5 , .08 , .2 , .92 , .7 , .13 , .2 , .5 , .7 , 1.05 , 1.0 , 1.0;
     .1 , .1 , .05 , .3 , .15 , .3 , .14 , .5 , .57 , .67 , .63 , .5 , .38 , .58];
C = [-1 , -4 ,  2 ,  3 , -8 ,   5 ,   6 ,   7;
     -2 ,  3 ,  4 , -6 ,  6 ,   7 ,   8 , -12;
     -3 , -5 ,  1 , -7 ,  2 ,   3 ,   4 , -13;
      3 ,  5 ,  6 ,  7 , -9 , -10 , -11 , -14];
H = [.22 , .12 , .3 , .2  , .07 , .3  , .2  , .3;
     .13 , .3  , .2 , .22 , .3  , .2  , .3  , .05;
     .05 , .15 , .2 , .16 , .2  , .2  , .2  , .12;
     .2  , .2  , .2 , .2  , .07 , .17 , .13 , .08];

figure();
plotGrid(C, H, G, B);

figure();
plotGridMesh(C, H, G, B);

%% 
Areal = 1e2 * ...
    [ 2.699300699300699                   0  -0.400000000000000                   0                   0                   0                   0                   0;
                      0   1.222222222222222  -0.158730158730159                   0  -0.285714285714286                   0                   0                   0;
     -0.250000000000000  -0.133333333333333   0.833333333333333  -0.200000000000000                   0  -0.250000000000000                   0                   0;
                      0                   0  -0.238095238095238   1.079545454545455                   0                   0  -0.277777777777778                   0;
                      0  -0.370370370370370                   0                   0   2.380952380952381  -0.180180180180180                   0                   0;
                      0                   0  -0.270270270270270                   0  -0.133333333333333   0.921568627450980  -0.200000000000000                   0;
                      0                   0                   0  -0.303030303030303                   0  -0.200000000000000   1.102564102564102  -0.133333333333333;
                      0                   0                   0                   0                   0                   0  -0.190476190476190   3.416666666666667];
Abreal = 1e2 * ...
    [-0.259740259740260  -0.439560439560440  -1.600000000000000                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0;
                      0                   0                   0  -0.396825396825397  -0.380952380952381                   0                   0                   0                   0                   0                   0                   0                   0                   0;
                      0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0;
                      0                   0                   0                   0                   0  -0.216450216450216  -0.347222222222222                   0                   0                   0                   0                   0                   0                   0;
                      0                   0                   0                   0                   0                   0                   0  -0.772200772200772  -1.058201058201058                   0                   0                   0                   0                   0;
                      0                   0                   0                   0                   0                   0                   0                   0                   0  -0.317965023847377                   0                   0                   0                   0;
                      0                   0                   0                   0                   0                   0                   0                   0                   0                   0  -0.466200466200466                   0                   0                   0;
                      0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0  -1.142857142857143  -0.833333333333333  -1.250000000000000];

[A, Ab] = shortleyWeller(C, H);
A = full(A);
Ab = full(Ab);
fA = max(abs((Areal - A))./A, [], 'all');
fAb = max(abs((Abreal - Ab)./Ab), [], 'all');

disp(['relativer Fehler in nichtnull Eintraegen von A ist: ' num2str(fA) '.']);
disp(['relativer Fehler in nichtnull Eintraegen von Ab ist: ' num2str(fAb) '.']);
