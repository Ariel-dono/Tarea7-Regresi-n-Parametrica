function Regresion_parametrica
    %N, number of observations
     Max = 1;
     N = 100;
     x = 0:Max/N:Max;
     
     %t_bold, set of target values
     [t, yReal] = getObservation(x);
     
     %M, degree of polynomial to fit
     M = 5;
	 %weights 
     w = getOptimumW(x, t, M);
     for i = 1:length(x)
        y_approx(i) = yF(x(i), w);
     end
     plotter(x,t,yReal,y_approx);
	 %measure the model error
     [Erms,error] = getError(yReal, y_approx);
     Erms
     lambda = 1/50;
     errorRegularizado = getErrorRegularizado(error,lambda,y_approx);
     errorRegularizado
     w = getOptimumRegularizadoW(x,t,M,lambda);
     %weights 
     for i = 1:length(x)
        y_approx(i) = yF(x(i), w);
     end
     plotter(x,t,yReal,y_approx);
     
end

function plotter(x,t,yReal,y_approx)
    figure;scatter(x, t);
    hold on;
    plot(x, y_approx);
    hold on;
    plot(x, yReal);
end
%RMS error function to evaluate
function [Erms,error] = getError(y, t)
    N = size(t);
    N = N(2);
    error = sum((y-t).^2)/2;
    Erms = sqrt((2*error)/N);
end

function[n] = distanciaEuclidiana(V)
 n = sqrt(sum(V.*V));
end

function result = getErrorRegularizado(error,lambda,w)
    result = error+(lambda/2)*distanciaEuclidiana(w');
end

function W = getOptimumRegularizadoW(x, t, M, lambda)
    N = size(x);
    N = N(2);
    X = zeros(N,M+1);
    for i=1:N
      for j=1:M+1
          X(i,j) = x(i)^(j-1);
      end
    end 
    W = pinv((X'*X)+(lambda*eye(M+1)))*X'*t';
end

%Least squares with first derivative equals zero
function W = getOptimumW(x, t, M)
   B = zeros(M+1,1);
   for i=1:M+1
       B(i,1) = sum(dot(x.^(i-1),t));
   end
   %B = sum(dot(x.^M,t));
   A = zeros(M+1,M+1);
   for i=1:M+1
      for j=1:M+1
          A(i,j) = sum(x.^((i-1)+(j-1)));
      end
   end
   A_1 = pinv(A);
   W = A_1*B;
end

%Phenomenom to analyze and adjust the model to
function [t, yReal] = getObservation(x)
    snr = 8;
    yReal = sin(2*pi*x);
    t = awgn(yReal, snr);
end

%Evaluaation of the polynomial model with weights w
function y = yF(x_val, w)
    y = 0;
    for i = 1:length(w)
        y = y + w(i) * (x_val.^(i - 1));
    end
end