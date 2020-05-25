
%%Ecuacion de la onda por diferencias finitas centradas
clc
clear all
h=input('Ingrese distancia de cada intervalo espacial:');
k=input('Ingrese distancia de cada intervalo temporal:');
a=input('Ingrese a=');
b=input('Ingrese b=');
T=input('Ingrese el tiempo final T=');
C=input('Ingrese el coeficiente C=');
%
p=input('Ingrese p punto de discontinuidad:');
fun1=input('Ingrese la función f definida en [a,p]:','s');
fun2=input('Ingrese la fución g definida en <p,b]:','s');
f1=inline(fun1);
f2=inline(fun2);
%
gun=input('Sea du/dt=g(x,0),ingrese: g(x)=','s');
g=inline(gun);
hun1=input('Sea u(a,t)=h1(t),ingrese h1(t)=','s');
hun2=input('Sea u(b,t)=h2(t),ingrese h2(t)=','s');
h1=inline(hun1);
h2=inline(hun2);
m=(b-a)/h;
n=T/k;
x=a+h:h:b-h;%columna
t=linspace(0,T,n+1);%fila
r=C^2*(k^2)/(h^2);
d=((p-a)/h);

u(1,1:n+1)=h1(t);
u(m+1,1:n+1)=h2(t);

u(2:d,1)=f1(a+h:h:p)';
u(d+1:m,1)=f2(p+h:h:b-h)';

u(2:d,2)=f1(a+h:h:p)'+k*g(a+h:h:p)';
u(d+1:m,2)=f2(p+h:h:b-h)'+k*g(p+h:h:b-h)';

A=2*(1-r)*eye(m-1)+r*diag(ones(1,m-2),+1)+r*diag(ones(1,m-2),-1);
for j=3:n+1
    H=[r*h1(t(j-1)),zeros(1,m-3),r*h2(t(j-1))];
    u(2:m,j)=A*u(2:m,j-1)-u(2:m,j-2)+H';
end
x=[a x b];
w=u';
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
   fprintf('Para t=%.4f tenemos:\n',T);
    for i=1:m+1
        fprintf('x(%3d)=%.5f:',i-1,x(i));
        fprintf('\t %25.10f',u(i,n+1));
        fprintf('\n');
    end
 %graficat
 Tfinal=T*ones(1,m+1);
 solfinal=w(n+1,:);
 plot3(x,Tfinal,solfinal);
 title('Solución aproximada para tiempo final');
 grid on;
