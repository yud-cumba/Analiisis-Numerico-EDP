%ECUACION DEL CALOR POR DIFERENCIAS FINITAS IMPLIC�TAS
clear all
clc
%******************************************************
fprintf('BIENVENIDOS AL PROGRAMA PARA DESARROLLAR\nLA ECUACION DEL CALOR POR DIFERENCIAS FINITAS IMPLIC�TAS\n');
disp('**Caso no homog�neo**');
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
fun=input('Ingrese la funci�n f definida en [a,b]:','s');
f=inline(fun);
%%%%%%%%%%%%%%%%%%%%
A=(1+2*r)*eye(m-1)+diag((-r)*ones(1,m-2),+1)+diag((-r)*ones(1,m-2),-1);
for i=2:m
    solfila(i-1)=f(x(i));
end

%solucion para t=0
h=input('Ingrese la funci�n en la frontera a:','s');
g=input('Ingrese la funci� en la frontera b:','s');
alfafrontera=inline(h);
betafrontera=inline(g);
sol=[alfafrontera(t(1)),solfila,betafrontera(t(1))];
 for j=2:n+1
     H=[-r*alfafrontera(t(j)),zeros(1,m-3),-r*betafrontera(t(j))];
     fila=inv(A)*(solfila'-H');
     solfila=fila';
     sol=[sol;alfafrontera(t(j)),solfila,betafrontera(t(j))];
 end
u=sol';
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
   fprintf('Para t=%.4f tenemos:\n',T);
    for i=1:m+1
        fprintf('x(%3d)=%.5f:',i-1,x(i));
        fprintf('\t %25.10f',u(i,n+1));
        fprintf('\n');
    end
 %graficat
 Tfinal=T*ones(1,m+1);
 solfinal=sol(n+1,:);
 plot3(x,Tfinal,solfinal);
 title('Soluci�n aproximada para tiempo final por el m�todo IMPL�CITO');
 grid on;
