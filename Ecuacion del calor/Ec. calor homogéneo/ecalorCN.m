%MÉTODO DE CRANK NICOLSON
clear all
clc
%******************************************************
fprintf('BIENVENIDOS AL PROGRAMA PARA DESARROLLAR\nLA ECUACION DEL CALOR POR CRANK NICOLSON\n');
fprintf('**Ecuación homogénea **\n');
w=input('Pulse una tecla para continuar...','s');
%******************************************************
h=input('Ingrese distancia de cada intervalo espacial:');
k=input('Ingrese distancia de cada intervalo temporal:');
a=input('Ingrese a=');
b=input('Ingrese b=');
T=input('Ingrese el tiempo final T=');
C=input('Ingrese el coeficiente C=');
m=(b-a)/h;
n=T/k;
x=linspace(a,b,m+1);
t=linspace(0,T,n+1);
r=C^2*k/(h^2);
ss=strcat('coeficiente r=',num2str(r));
disp(ss); 
fun=input('Ingrese la función f:','s');
f=inline(fun);

%%%%%%%%%%%%%%%%%%%%
A=(1+r)*eye(m-1)+diag((-r/2)*ones(1,m-2),+1)+diag((-r/2)*ones(1,m-2),-1);
B=(1-r)*eye(m-1)+diag((r/2)*ones(1,m-2),+1)+diag((r/2)*ones(1,m-2),-1);

for i=2:m
    solfila(i-1)=f(x(i));
end
%solucion para t=0
sol=[f(x(1)),solfila,f(x(m+1))];
 for j=2:n+1
     fila=inv(A)*B*solfila';
     solfila=fila';
     sol=[sol;0,solfila,0];
 end
u=sol';
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('Para t=%.4f tenemos:\n',t(n+1));
    for i=1:m+1
        fprintf('x(%3d)=%.5f:',i-1,x(i));
        fprintf('\t %25.10f',u(i,n+1));
        fprintf('\n');
    end
 %graficat
 Tfinal=T*ones(1,m+1);
 solfinal=sol(n+1,:);
 plot3(x,Tfinal,solfinal);
 title('Solución aproximada para tiempo final por el método CRANK NICOLSON');
 grid on;