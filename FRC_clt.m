clc
q1=10^9.*[181.8,2.897,0;2.897,10.35,0;0,0,7.17];                              %Glass/epoxy stiffness matrix
q2=10^3.*[137.69,3.34,0;3.34,11.93,0;0,0,6.6];                                               %Carbon/epoxy stiffness matrix
n=input('enter no of layers:');                      %number of layers of lalminate
t=input('enter the thickness of laminate:');         %thickness of each layer
h=t/n;                                               %calculating the height of each layer
Nx=input('enter Nx:');                                  %enter loads and moments
Ny=input('enter Ny:');
Nxy=input('enter Nxy:');
Mx=input('enter Mx:');
My=input('enter My:');
Mxy=input('enter Mxy:');
N=[Nx;Ny;Nxy;Mx;My;Mxy];
syms sigma
total_matrix=zeros(6);    
z=[];
angle_total=[];                                                            
for i=1:n
  angle=input('enter the angle of layer:');         %enter the angle of that layer in degree
  p=input('enter the material of the layer: ');     %determining the material of the composite by entering 1or 2
  if p==1
      q=q1;
  else
      q=q2;
  end    
  a=(angle*3.14)/180;                                  %convert angle into radians
  angle_total(i)=a;                                      
  h1=-t/2+(i-1)*h;                                     %determine the z-axis values corresponding to layer
  h2=-t/2+i*h;
  if i==1
      z(i)=h1;
  end
  z(i+1)=h2;                                            %storing z values in a matrix
  q_actual=Qmatrix_calculation(q,a)                    %calling a function to calculate Qmatrix for given laminate 
  matrix=matrix_calculation(h1,h2,q_actual);            %calling a function to calculate ABD matrix of each layer
  total_matrix=total_matrix+matrix;                     %summation of all ABD matrix of each layer
end  
disp(z)                                                 
disp(angle_total)
disp(total_matrix)
g=inv(total_matrix);                                    %inversion of ABD matrix
disp(g);
strain_matrix=g*N;                                       %calculation of strain matrix
disp(strain_matrix)
e=[strain_matrix(1);strain_matrix(2);strain_matrix(3)];  
k=[strain_matrix(4);strain_matrix(5);strain_matrix(6)];
stress=[];
for i=1:n
    q_actual=Qmatrix_calculation(q,angle_total(i));        
    for j=i:i+1
        h=e+(z(j)*k);
        stress(:,j)=q_actual*h;                          %stress calculation from the Q and strains of laminates for each value of z
    end
end
disp(stress)
plot(z,stress(1,:),'k')
hold on
plot(z,stress(2,:),'b')
hold on 
plot (z,stress(3,:),'r')
xlabel('Thickness(t)mm')
ylabel('Stress(Mpa)')
legend('sigmax','sigmay','sigmaxy')
      


function q_actual=Qmatrix_calculation(q,theeta)
  m=cos(theeta);
  n=sin(theeta);
  R=[1,0,0;0,1,0;0,0,2];                                 
  T=[m^2,n^2,2*m*n;n^2,m^2,-2*n*m;-n*m,m*n,m^2-n^2];
  inv_R=inv(R);
  inv_T=inv(T);
  q_actual=inv_T*q*R*T*inv_R;                               %transposed Qmatrix for given lamina material and angle
end  

function matrix=matrix_calculation(h1,h2,q_actual)
  
  A=zeros(3);
  B=zeros(3);
  D=zeros(3);
  
  
  for i=1:3
    for j=1:3
      A(i,j)=A(i,j)+q_actual(i,j)*(h2-h1);                %A matrix determination
    end  
  end
  for i=1:3
    for j=1:3
      B(i,j)=B(i,j)+(0.5*q_actual(i,j)*(h2^2-h1^2));      %B matrix determination
    end
  end    

  for i=1:3
    for j=1:3
      D(i,j)=D(i,j)+(0.33*q_actual(i,j)*(h2^3-h1^3));   %D matrix determination
    end
  end  
  matrix=[A,B;B,D];
  disp(matrix)                                 
end  


