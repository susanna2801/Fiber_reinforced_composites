
clc
q1=10^3.*[47.709,10.49,0;10.49,20.386,0;0,0,6];                              %Glass/epoxy stiffness matrix
q2=10^3.*[137.69,3.34,0;3.34,11.93,0;0,0,6.6];                                               %Carbon/epoxy stiffness matrix

n=input('enter no of layers:');                       %number of layers of lalminate    
t=input('enter the thickness of laminate:');          %thickness of each layer
h=t/n;                                                %calculating the height of each layer
alpha1=10^-6.*[5.41;9.92;9.92];                                        %coefficient of thermal expansion for material
alpha2=10^-6.*[0.9;27;27];
Tref=input('Enter Tref:');                            % define the reference temperature and operational temperature
T=input('Enter Toper:');
tem=T-Tref;
N=zeros(6,1);

total_matrix=zeros(6);                            %creating a matrix of zeros to add ABD matrix
z=[];
angle_total=[];
for i=1:n                                            %using for loop to calculate ABD matrix for each layer
  angle=input('enter the angle of layer:');         %enter the angle of that layer
  a=(angle*3.14)/180;
  p=input('enter the material of the layer: ');     %determining the material of the composite by entering 1or 2
  if p==1
      q=q1;
      alpha=alpha1;
  else
      q=q2;
      alpha=alpha2;
  end    
  angle_total(i)=a;                                   %convert angle into radians
  h1=-t/2+(i-1)*h;                                    %determine the z-axis values corresponding to layer
  h2=-t/2+i*h;
  if i==1
      z(i)=h1;
  end
  z(i+1)=h2;                                        %storing z values in a matrix
  alpha_actual=alpha_calculation(alpha,a);
  thermal=thermal_loads(h1,h2,alpha_actual,T,Tref);        %calling a function to calculate thermal loads of each layer     
  q_actual=Qmatrix_calculation(q,a);                %calling a function to calculate Qmatrix for given laminate
  matrix=matrix_calculation(h1,h2,q_actual);         %calling a function to calculate ABD matrix of each layer
  total_matrix=total_matrix+matrix;                 %summation of all ABD matrix of each layer
  N=N+thermal;                                     %storing thermal loads in N matrix  
end  
disp(z)
disp(angle_total);
disp(total_matrix);
g=inv(total_matrix);                                  
%inversion of ABD matrix
strain_matrix=g*N;                     %calculation of strain matrix
disp(strain_matrix);      
e=[strain_matrix(1);strain_matrix(2);strain_matrix(3)];
k=[strain_matrix(4);strain_matrix(5);strain_matrix(6)];
stress=[];
for i=1:n
    l=input('Enter the material:');             %enter either 1 or 2 for glass epoxy or carbon epoxy for stress calculation
    if l==1
        Q=q1;
    else
        Q=q2;
    end    
    alpha_actual=alpha_calculation(alpha,angle_total(i));
    ethermal=[alpha_actual(1)*t;alpha_actual(1)*t;alpha_actual(1)*t];
    q_actual=Qmatrix_calculation(Q,angle_total(i));     
    for j=i:i+1
        h=e+(z(j)*k)-ethermal;               %summing of strain and curvature
        stress(:,j)=q_actual*h;               %stress calculation from strain and Qmatrix
    end
end

disp(stress)
plot(z,stress(1,:))
hold on
plot(z,stress(2,:))
hold on
plot(z,stress(3,:))
xlabel('Thickness(t)mm');
ylabel('Stress(Mpa)');
legend('sigmax','sigmay','sigmaxy');

function alpha_actual=alpha_calculation(alpha,theeta)
  m=cos(theeta);
  n=sin(theeta);
  alpha_actual(1)=alpha(1)*m^2+alpha(2)*n^2;
  alpha_actual(2)=alpha(1)*n^2+alpha(2)*m^2;                    %Transformation of coefficient of linear expansion
  alpha_actual(3)=alpha(2)*n*m-alpha(1)*m*n;

end  

      


function q_actual=Qmatrix_calculation(q,theeta)
  m=cos(theeta);
  n=sin(theeta);
  R=[1,0,0;0,1,0;0,0,2];
  T=[m^2,n^2,2*m*n;n^2,m^2,-2*n*m;-n*m,m*n,m^2-n^2];
  inv_R=inv(R);
  inv_T=inv(T);
  q_actual=inv_T*q*R*T*inv_R;                    %transposed Qmatrix for given lamina material and angle
end  

function matrix=matrix_calculation(h1,h2,q_actual)
  
  A=zeros(3);
  B=zeros(3);
  D=zeros(3);
  
  
  for i=1:3
    for j=1:3
      A(i,j)=A(i,j)+q_actual(i,j)*(h2-h1);          %A matrix determination
    end  
  end
  for i=1:3
    for j=1:3
      B(i,j)=B(i,j)+(0.5*q_actual(i,j)*(h2^2-h1^2));      %B matrix determination
    end
  end    

  for i=1:3
    for j=1:3
      D(i,j)=D(i,j)+(0.33*q_actual(i,j)*(h2^3-h1^3));   %C matrix determination
    end
  end  
  matrix=[A,B;B,D];
  disp(matrix)
end 

function thermal=thermal_loads(h1,h2,alpha,T,Tref)
  Nx=alpha(1)*(h2-h1)*(T-Tref);
  Ny=alpha(2)*(h2-h1)*(T-Tref);
  Nxy=alpha(3)*(h2-h1)*(T-Tref);
  Mx=0.5*alpha(1)*(h2^2-h1^2)*(T-Tref);                %calculating thermal loads from alpha,Tref,Toper and thickness of each layer
  My=0.5*alpha(2)*(h2^2-h1^2)*(T-Tref);
  Mxy=0.5*alpha(3)*(h2^2-h1^2)*(T-Tref);
  thermal=[Nx;Ny;Nxy;Mx;My;Mxy];
end  

