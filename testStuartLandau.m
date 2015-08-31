% Try DMD/cDMD with MW's basis functions
clear all
close all
dt = 0.01
t= 0:dt:10;
r0 = 2;

mu = 1;
th0 = 0;
gamma = 1;
beta = 0.;


r2 = r0*(mu./(r0^2+(mu-r0^2)*exp(-2*mu*t))).^0.5;
th = th0 + (gamma-mu*beta)*t+beta*log(r2/r0);

minr = -1;
maxr = 0;

minn = -3;
maxn = 3;
clear i
clear Y
Nrows = (maxn-minn+1)*(maxr-minr+1);
Y = zeros(Nrows,length(t));
%for tt = 1:length(t)
    k= 1;
    for rr = minr:maxr
        for nn = minn:maxn
            Y(k,:) = r2.^(2*rr).*exp(1i*nn*th);
            k=k+1;
        end
    end
%end
order = Nrows;
% do for a range of \mu

Mus = 0.1:0.1:0.5;
Betas = 0:0.01:0.2;
Eigs = zeros(order,length(Mus));
EigsSort = zeros(order,length(Mus));
for ii = 1:length(Betas)
    %mu = Mus(ii)
    beta = Betas(ii)
  %r0 = Mus(ii)
    r2 = r0*(mu./(r0^2+(mu-r0^2)*exp(-2*mu*t))).^0.5;
    th = th0 + (gamma-mu*beta)*t+beta*log(r2/r0);

  %  for tt = 1:length(t)
      %  tt
      k=1;
        for rr = minr:maxr
            for nn = minn:maxn
                Y(k,:) = r2.^(2*rr).*exp(1i*nn*th);
                k=k+1;
            end
        end
   % end
    
    [PsiE, LambdaE, U, S, V,] = DMDext(Y,order);
    Eigs(:,ii) = diag(LambdaE);
    [~,EigsSortRealInd] = sort(real(Eigs(:,ii)));
    EigsSortReal = Eigs(EigsSortRealInd,ii);
    [EigsSortImag1,EigsSortImag1Ind]  = sort(imag(EigsSortReal(1:end/2)));
    [EigsSortImag2,EigsSortImag2Ind]  = sort(imag(EigsSortReal((end/2+1):end)));
    EigsSort(:,ii) = [EigsSortReal(EigsSortImag1Ind);EigsSortReal(EigsSortImag2Ind+(maxn-minn+1))];
    
end



figure
eignumber = 1;
for ii = 1:length(Betas)
    plot(real(log((Eigs(:,ii)))/dt),imag(log((Eigs(:,ii)))/dt),'x','color',[ii/length(Betas),0,ii/length(Betas)])
    hold on
end
grid on

%%
figure
%eignumber = 8
for eignumber = 1:14
%for ii = 1:length(Betas)
    plot(real(log((EigsSort(eignumber,:)))/dt),imag(log((EigsSort(eignumber,:)))/dt),'x')
    hold on
%end
end
grid on

