%%This fuction tests my analysis of the the elipsoidal transformations and
%%derivatives. We attempt to confirm how the representation of a polyhedra
%%changes under the proposed elipsoidal tranforms 
close all

r = 1.5;
sensor_angle = pi/5;
triangle = [0 cos(sensor_angle) cos(sensor_angle) 0;
            0 -sin(sensor_angle) sin(sensor_angle) 0;]*r;

interior_point = [.01; 0];
up_point = [.01; 1];
down_point = [.01; -1];
out_point = [r+.01; 0];
    
%%One or both of the normals here may be outward facing,
%%need to check this 
normals = [-sin(sensor_angle) -sin(sensor_angle)  1;
           -cos(sensor_angle) cos(sensor_angle)   0;];
        
plot_normals = [ 0 sin(sensor_angle) 0 -sin(sensor_angle) 0 -1;
            0 cos(sensor_angle) 0 cos(sensor_angle)  0  0;];
        
b = [0;0;cos(sensor_angle)*r;];
        
rot_ang = pi/3;         
        
Q = [1 0;
      0 2;];
  
Rtest = [cos(.3), -sin(.3);
      sin(.3),  cos(.3);];
  
Q = Rtest*Q;
  
Qinv = inv(Q);


%%The identity chol(Q) = inv(chol(Qinv))' is only true for an axis-aligned
%%elipse. This means that we may have to rotate not by the angle between 
%%poses, but rather the required angle to make Q axis aligned  
chol(Q')
inv(chol(Qinv))'
        
Rk = [cos(rot_ang), -sin(rot_ang);
      sin(rot_ang),  cos(rot_ang);];

rot_tri = Rk*triangle;  
rot_plot_normals = Rk*plot_normals;
rot_normals = Rk*normals;
rot_point = Rk*interior_point;
rot_up = Rk*up_point;
rot_down = Rk*down_point;
rot_out = Rk*out_point;


scew_tri = chol(Qinv)*rot_tri;
plot_scew_norm = inv(chol(Qinv))'*rot_plot_normals;
scew_norm = inv(chol(Qinv))'*rot_normals;
scew_point = inv(chol(Qinv))'*rot_point;
scew_up = Rk*rot_up;
scew_down = Rk*rot_down;
scew_out = Rk*rot_out;

b_new = b;
for i =1:length(b)
b_new(i) = b(i)/(norm(scew_norm(:,i)));

end

normals'*interior_point
normals'*up_point
normals'*down_point
normals'*out_point
b

scew_norm'*scew_point
scew_norm'*scew_up
scew_norm'*scew_down
scew_norm'*scew_out
b_new


plot(rot_tri(1,:), rot_tri(2,:));
hold on 
plot(-rot_tri(1,:), -rot_tri(2,:));

plot(rot_plot_normals(1,:), rot_plot_normals(2,:), 'g');

hold off
axis([-2 2 -2 2])


figure 

plot(scew_tri(1,:), scew_tri(2,:));
hold on 
plot(-scew_tri(1,:), -scew_tri(2,:));

plot(plot_scew_norm(1,:), plot_scew_norm(2,:), 'g');

hold off
axis([-2 2 -2 2])

plot_scew_norm(:,2)'*scew_tri(:,2)
plot_scew_norm(:,4)'*scew_tri(:,3)