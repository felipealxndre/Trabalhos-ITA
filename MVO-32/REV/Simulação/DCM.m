function C = DCM(i_dir, a_rad) %matriz de cossenos diretores

%3:rotacao em z; 2:rotacao em y; 1:rotacao em x

switch i_dir
    case 1
        C = [1 0 0
             0 cos(a_rad) sin(a_rad)
             0 -sin(a_rad) cos(a_rad)
             ];

    case 2
        C = [cos(a_rad) 0 -sin(a_rad)
            0 1 0
            sin(a_rad) 0 cos(a_rad)
             ];
    
    case 3
           C = [cos(a_rad) sin(a_rad) 0
               -sin(a_rad) cos(a_rad) 0
               0 0 1
               ];
    otherwise
       error('i_dir must be 1, 2 or 3!')
end