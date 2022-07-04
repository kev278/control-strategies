% Code by Lauren Stanley
% Group Project 

function [x_desired,y_desired,z_desired]=TrajectoryGeneratorv0()
    plotOn = false;
    syms t;    
    startvelocity = 0;
    endvelocity = 0;
    startaccel = 0;
    endaccel = 0;
    p0 = [0,0,0];
    p1 = [0,0,1];
    p2 = [1,0,1];
    p3 = [1,1,1];
    p4 = [0,1,1];
    p5 = [0,0,1];
    t0 = 0;
    t1 = 5;
    t2 = 20;
    t3 = 35;
    t4 = 50;
    t5 = 65;
    
    x01 = quinticpoly(t0,t1,p0(1),p1(1),startvelocity,endvelocity,startaccel,endaccel);
    y01 = quinticpoly(t0,t1,p0(2),p1(2),startvelocity,endvelocity,startaccel,endaccel);
    z01 = quinticpoly(t0,t1,p0(3),p1(3),startvelocity,endvelocity,startaccel,endaccel);
    
    x12 = quinticpoly(t1,t2,p1(1),p2(1),startvelocity,endvelocity,startaccel,endaccel);
    y12 = quinticpoly(t1,t2,p1(2),p2(2),startvelocity,endvelocity,startaccel,endaccel);
    z12 = quinticpoly(t1,t2,p1(3),p2(3),startvelocity,endvelocity,startaccel,endaccel);
    
    x23 = quinticpoly(t2,t3,p2(1),p3(1),startvelocity,endvelocity,startaccel,endaccel);
    y23 = quinticpoly(t2,t3,p2(2),p3(2),startvelocity,endvelocity,startaccel,endaccel);
    z23 = quinticpoly(t2,t3,p2(3),p3(3),startvelocity,endvelocity,startaccel,endaccel); 
    
    x34 = quinticpoly(t3,t4,p3(1),p4(1),startvelocity,endvelocity,startaccel,endaccel);
    y34 = quinticpoly(t3,t4,p3(2),p4(2),startvelocity,endvelocity,startaccel,endaccel);
    z34 = quinticpoly(t3,t4,p3(3),p4(3),startvelocity,endvelocity,startaccel,endaccel);
    
    x45 = quinticpoly(t4,t5,p4(1),p5(1),startvelocity,endvelocity,startaccel,endaccel);
    y45 = quinticpoly(t4,t5,p4(2),p5(2),startvelocity,endvelocity,startaccel,endaccel);
    z45 = quinticpoly(t4,t5,p4(3),p5(3),startvelocity,endvelocity,startaccel,endaccel);
    
    x_desired = piecewise(0<t<=t1,x01,t1<t<=t2,x12,t2<t<=t3,x23,t3<t<=t4,x34,t4<t<=t5,x45);
    y_desired = piecewise(0<t<=t1,y01,t1<t<=t2,y12,t2<t<=t3,y23,t3<t<=t4,y34,t4<t<=t5,y45);
    z_desired = piecewise(0<t<=t1,z01,t1<t<=t2,z12,t2<t<=t3,z23,t3<t<=t4,z34,t4<t<=t5,z45);
    
    if plotOn
        clf;
        
        subplot(3,1,1)
        fplot(x_desired,[t0,t5],'LineWidth',2.0);
        xlabel('time [s]'), ylabel('position [m]'), title('X Position');
        
        subplot(3,1,2)
        fplot(y_desired,[t0,t5],'LineWidth',2.0);
        xlabel('time [s]'), ylabel('position [m]'), title('Y Position');
        
        subplot(3,1,3)
        fplot(z_desired,[t0,t5],'LineWidth',2.0);
        xlabel('time [s]'), ylabel('position [m]'), title('Z Position');
        
        pause;
        clf;
        
        subplot(3,1,1)
        fplot(x_dot_desired,[t0,t5],'LineWidth',2.0);
        xlabel(['time [s]']), ylabel('velocity [m/s]'), title('X Velocity');
        
        subplot(3,1,2)
        fplot(y_dot_desired,[t0,t5],'LineWidth',2.0);
        xlabel('time [s]'), ylabel('velocity [m/s]'), title('Y Velocity');
        
        subplot(3,1,3)
        fplot(z_dot_desired,[t0,t5],'LineWidth',2.0);
        xlabel('time [s]'), ylabel('velocity [m/s]'), title('Z Velocity');
        
        pause;
        clf;
        
        subplot(3,1,1)
        fplot(x_double_dot_desired,[t0,t5],'LineWidth',2.0);
        xlabel('time [s]'), ylabel('acceleration [m/s^2]'), title('X Acceleration');
        
        subplot(3,1,2)
        fplot(y_double_dot_desired,[t0,t5],'LineWidth',2.0);
        xlabel('time [s]'), ylabel('acceleration [m/s^2]'), title('Y Acceleration');
        
        subplot(3,1,3)
        fplot(z_double_dot_desired,[t0,t5],'LineWidth',2.0);
        xlabel('time [s]'), ylabel('acceleration [m/s^2]'), title('Z Acceleration');
    end
end

function q = quinticpoly(t0,tf,q0,qf,qd0,qdf,qdd0,qddf)
    A1 = [1 t0 t0^2 t0^3 t0^4 t0^5;...
    0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;...
    0 0 2 6*t0 12*t0^2 20*t0^3;...
    1 tf tf^2 tf^3 tf^4 tf^5;...
    0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;...
    0 0 2 6*tf 12*tf^2 20*tf^3];
    b1 = [q0 qd0 qdd0 qf qdf qddf]';
    
    A = A1\b1;
    syms t;
    q(t) = A(1) + A(2) * t + A(3) * t^2 + A(4) * t^3 + A(5) * t^4 + A(6) * t^5;
end
