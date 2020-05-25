%Ecuacion de la onda por diferencias finitas centradas
clc
clear all
m=input('Ingrese numero intervalos variable espacial:');
n=input('Ingrese numero intervalos variable temporal:');
a=input('Ingrese a=');
b=input('Ingrese b=');
T=input('Ingrese el tiempo final T=');
C=input('Ingrese el coeficiente C=');
fun=input('Sea u(x,0)=f(x),ingrese: f(x)=','s');
gun=input('Sea du/dt=g(x,0),ingrese: g(x)=','s');
f=inline(fun);
g=inline(gun);
u=zeros(m+1,n+1);
h=(b-a)/m;
k=T/n;
r=C^2*(k^2)/(h^2);
x=a+h:h:b-h;% x\in <a,b>
t=linspace(0,T,n+1);
u(2:m,1)=f(x);
u(2:m,2)=(1-r)*f(x)'+k*g(x)'+r/2*(u(1:m-1,1)+u(3:m+1,1));

A=2*(1-r)*eye(m-1)+r*diag(ones(1,m-2),+1)+r*diag(ones(1,m-2),-1);
for j=3:n+1
    u(2:m,j)=A*u(2:m,j-1)-u(2:m,j-2);
end
x=[a x b];
w=u';
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    for i=1:m+1
        fprintf('x(%3d)=%.4f',i-1,x(i));
        for j=1:n+1
        fprintf('\t %25.10f',u(i,j));
        end
        fprintf('\n');
    end
surf(x,t,w);
xlabel('x');
ylabel('t');