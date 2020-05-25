
clear all
clc
%******************************************************
fprintf('BIENVENIDOS AL PROGRAMA PARA DESARROLLAR\nLA ECUACION DEL CALOR POR CRANK NICOLSON\n');
fprintf('**Ecuaci�n no homog�nea con f(x)\n que tiene una regla de correspondencia por tramos**\n');
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
fun=input('Ingrese la funci�n f:','s');
f=inline(fun);
%%%%%%%%%%%%%%%%%%%%%
p=input('Ingrese el punto de discontinuidad p:');
fun1=input('Ingrese la funci�n f definida en [a,p>:','s');
fun2=input('Ingrese la fuci�n g definida en [p,b]:','s');
f1=inline(fun1);
f2=inline(fun2);

%%%%%%%%%%%%%%%%%%%%
A=(1+r)*eye(m-1)+diag((-r/2)*ones(1,m-2),+1)+diag((-r/2)*ones(1,m-2),-1);
B=(1-r)*eye(m-1)+diag((r/2)*ones(1,m-2),+1)+diag((r/2)*ones(1,m-2),-1);


d=((p-a)/h);
%+1 si [a,p]
for i=2:d
    solfila(i-1)=f1(x(i));
end
for i=d+1:m
    solfila(i-1)=f2(x(i));
end

%solucion para t=0
sol=[f(x(1)),solfila,f(x(m+1))];
h=input('Ingrese la funci�n en la frontera a:','s');
g=input('Ingrese la funci� en la frontera b:','s');
alfafrontera=inline(h);
betafrontera=inline(g);
sol=[alfafrontera(t(1)),solfila,betafrontera(t(1))];
 for j=2:n+1
     H1=[-r/2*alfafrontera(t(j-1)),zeros(1,m-3),-r/2*betafrontera(t(j-1))];
     H2=[r/2*alfafrontera(t(j)),zeros(1,m-3),r/2*betafrontera(t(j))];
     fila=inv(A)*(B*solfila'+H2'-H1');
     solfila=fila';
     sol=[sol;alfafrontera(t(j)),solfila,betafrontera(t(j))]; 
 end
sol;
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
 title('Soluci�n aproximada para tiempo final por el m�todo CRANK NICOLSON');
 grid on;