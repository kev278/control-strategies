% Code by Lauren Stanley
% Group Project

function [x_desired,y_desired,z_desired,x_dot_desired,y_dot_desired,z_dot_desired,x_double_dot_desired,y_double_dot_desired,z_double_dot_desired]=TrajectoryGenerator(t)  
    t0 = 0;
    t1 = 5;
    t2 = 20;
    t3 = 35;
    t4 = 50;
    t5 = 65;    
    if t0<t&&t<=t1
        x_desired = 0;
        y_desired = 0;
        z_desired = (6*t^5)/3125 - (3*t^4)/125 + (2*t^3)/25 - (5404319552844595*t^2)/649037107316853453566312041152512 - t/649037107316853453566312041152512;
        x_dot_desired = 0;
        y_dot_desired = 0;
        z_dot_desired = (6*t^4)/625 - (12*t^3)/125 + (6*t^2)/25 - (5404319552844595*t)/324518553658426726783156020576256 - 1/649037107316853453566312041152512;
        x_double_dot_desired = 0;
        y_double_dot_desired = 0;
        z_double_dot_desired = (24*t^3)/625 - (36*t^2)/125 + (12*t)/25 - 5404319552844595/324518553658426726783156020576256;
    elseif t1<t&&t<=t2        
        x_desired = (4664065662093481*t^5)/590295810358705651712 - t^4/2025 + (22*t^3)/2025 - (8*t^2)/81 + (32*t)/81 - 47/81;
        y_desired = 0;
        z_desired = 1;
        x_dot_desired = (23320328310467405*t^4)/590295810358705651712 - (4*t^3)/2025 + (22*t^2)/675 - (16*t)/81 + 32/81;
        y_dot_desired = 0;
        z_dot_desired = 0;
        x_double_dot_desired = (23320328310467405*t^3)/147573952589676412928 - (4*t^2)/675 + (44*t)/675 - 16/81;
        y_double_dot_desired = 0;
        z_double_dot_desired = 0;
    elseif t2<t&&t<=t3    
        x_desired = 1;
        y_desired = (4664065662093517*t^5)/590295810358705651712 - (11*t^4)/10125 + (118*t^3)/2025 - (616*t^2)/405 + (1568*t)/81 - 7808/81;
        z_desired = 1;
        x_dot_desired = 0;
        y_dot_desired = (23320328310467585*t^4)/590295810358705651712 - (44*t^3)/10125 + (118*t^2)/675 - (1232*t)/405 + 1568/81;
        z_dot_desired = 0;
        x_double_dot_desired = 0;
        y_double_dot_desired = (23320328310467585*t^3)/147573952589676412928 - (44*t^2)/3375 + (236*t)/675 - 1232/405;
        z_double_dot_desired = 0;
    elseif t3<t&&t<=t4
        x_desired = - (4664065662091993*t^5)/590295810358705651712 + (7743077759332409*t^4)/4611686018427387904 - (5088511578973043*t^3)/36028797018963968 + (413524965784659*t^2)/70368744177664 - (8513749295566505*t)/70368744177664 + 8687499281190311/8796093022208;
        y_desired = 1;
        z_desired = 1;
        x_dot_desired = - (23320328310459965*t^4)/590295810358705651712 + (7743077759332409*t^3)/1152921504606846976 - (15265534736919129*t^2)/36028797018963968 + (413524965784659*t)/35184372088832 - 8513749295566505/70368744177664;
        y_dot_desired = 0;
        z_dot_desired = 0;
        x_double_dot_desired = - (23320328310459965*t^3)/147573952589676412928 + (23229233277997227*t^2)/1152921504606846976 - (15265534736919129*t)/18014398509481984 + 413524965784659/35184372088832;
        y_double_dot_desired = 0;
        z_double_dot_desired = 0;
    elseif t4<t&&t<=t5
        x_desired = 0;
        y_desired = - (4664065662091041*t^5)/590295810358705651712 + (1309491091651539*t^4)/576460752303423488 - (2339647806415457*t^3)/9007199254740992 + (4156099656120621*t^2)/281474976710656 - (7340936892604373*t)/17592186044416 + 2579651104916761/549755813888;
        z_desired = 1;
        x_dot_desired = 0;
        y_dot_desired = - (23320328310455205*t^4)/590295810358705651712 + (1309491091651539*t^3)/144115188075855872 - (7018943419246371*t^2)/9007199254740992 + (4156099656120621*t)/140737488355328 - 7340936892604373/17592186044416;
        z_dot_desired = 0;
        x_double_dot_desired = 0;
        y_double_dot_desired = - (23320328310455205*t^3)/147573952589676412928 + (3928473274954617*t^2)/144115188075855872 - (7018943419246371*t)/4503599627370496 + 4156099656120621/140737488355328;
        z_double_dot_desired = 0;
    else
        x_desired = 0;
        y_desired = 0;
        z_desired = 0;
        x_dot_desired = 0;
        y_dot_desired = 0;
        z_dot_desired = 0;
        x_double_dot_desired = 0;
        y_double_dot_desired = 0;
        z_double_dot_desired = 0;
    end
end
