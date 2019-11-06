% Targeted importance sampling-single importance sampling
% Sum of gaussians
clear 
% close all

tic
Ovar=1; % variance of original distribution

gaussian=@(x,mu,v) 1/sqrt(2*pi*v).*exp(-(x-mu).^2/(2*v)); 
N=10; % number of gaussians in sum
m=10;%-15:1:15; % biasing target
ev=[N,N/2,2,1];
J=length(ev);
mu=m/N;

numSamples=50000; 
seeds=1:numSamples*J; % for replicability
numBins=101;
M=numSamples;
% locs are bin centers
locs=linspace(-25,25,numBins);
binwidth = locs(2)-locs(1);
% Phat=zeros(size(locs));
Phat=zeros(J,length(locs)); 
Shat=zeros(J,length(locs)); 

% for coefficient of variation
bintracker=zeros(size(locs));
binmeans=zeros(size(locs));

%% no off

for j=1:J 
    bin=zeros(M,1);
    binweight=bin;

    % set up so that rng goes from 1:numSamples
%     if j>1
%         base=sum(M(1:j-1));
%     else
%         base=0;
%     end
    base=(j-1)*M;
    for k=1:M
        allmu=mu*ones(N,1);
        % draw sample
        rng(seeds(base+k))
        x=sqrt(Ovar)*randn(N,1); z=0;
        for nn=1:N
            if mod(nn,ev(j))==0 && (nn<N || ev(j)==1)
                newmu=(m-z)/(N-nn+1);
                allmu(nn:end)=newmu*ones(size(allmu(nn:end)));
            end
            z=z+x(nn)+allmu(nn);
        end
        xstar=x+allmu;
%         z=sum(xstar);
        
        % binning
        finder=abs(locs-z);
        bin(k)=find(finder==min(finder),1); %targeted
        % likelihood ratio-each x is iid gaussian
%         biased=zeros(1,J); 
%         for ll=1:J
%            biased(ll)=prod(gaussian(xstar,mu(ll),Ovar),1);
%         end
        biased=prod(gaussian(xstar,allmu,Ovar),1);
        invlike=biased./prod(gaussian(xstar,0,Ovar)); 
        binweight(k)=1/sum(M'.*invlike);
%         likes(k,:)=1./invlike;
    end
%     nooff(j).z=z;
%     nooff(j).likes=likes;
    for k=1:M
        % Now do the importance sampling sum
        Phat(j,bin(k))=Phat(j,bin(k))+binweight(k);%*1/sum(invlike);  
    end
    Phat(j,:)=Phat(j,:)./binwidth;
    for k=1:M
        thisbin=bin(k);
        temp=Phat(j,:);
        mybinweight=M*binweight(k);
        temp(thisbin)=Phat(j,thisbin)-mybinweight;
        Shat(j,:)=Shat(j,:)+temp.^2;
    end
%     nooff(j).bin=bin;
end

% importancesampledpdfintegral=1/binwidth*trapz(locs,Sample)
% exactpdfintegral=trapz(locs,exact(locs,10*myvar))
% mean=1/(max(locs)-min(locs))*trapz(locs,Sample)
% Phatn=sum(Phat,1); 
% varn=sum(1./(M.*(M-1)).*Shat,1); 
varn=1./(M.*(M-1)).*Shat;
covn=sqrt(varn)./Phat; 

%%
% binmc = zeros(1,numSamples);
% zz = zeros(1,numSamples);
% % Straight monte-carlo
% for i=1:numSamples
%     xmc = sqrt(Ovar)*randn(N,1);
%     zz(i) = sum(xmc);
%     [~,binmc(i)] = min(abs(locs - zz(i)));
% end
% 
% mcP = zeros(size(Phatn)); mcV = mcP;
% for j = 1:numBins
%     inds = find(binmc == j);
%     mcP(j) = length(inds)/numSamples/binwidth;
%     mcV(j) = sum(((binmc==j)-mcP(j)).^2);
% end
% 
% mcV = 1/(numSamples*(numSamples-1))*mcV;
% 
% mcCV = sqrt(mcV)./mcP;

%% Plotting
figure(1)
subplot(2,1,1)
semilogy(locs,gaussian(locs,0,N*Ovar),'-g','linewidth',2)
legendinfo{1}='Exact';
hold on
linestyle={'-','--','-.',':'};
semilogy(locs,Phat(1,:),linestyle{1},'linewidth',2)
legendinfo{1+1}='Regular ISMC';
for h=2:J-1
    semilogy(locs,Phat(h,:),linestyle{h},'linewidth',2)
    legendinfo{h+1}=sprintf('Recalculated every %d steps',ev(h));
end
semilogy(locs,Phat(J,:),linestyle{J},'linewidth',2)
legendinfo{J+1}='Recalculated every step';
plot(m*ones(100,1),linspace(1e-100,1),'--','color',[192,192,192]/255,'linewidth',2)
title([num2str(numSamples),' samples of Gaussian random walk'])
ylabel('Probability')
legend(legendinfo,'location','southwest')
axis([-5 25 1e-14 1])
hold off%,['target ',num2str(m(1))],['target ',num2str(m(2))],...
%     ['target ',num2str(m(3))])

% figure(10)
subplot(2,1,2)
% semilogy(locs,mcCV,'--b',locs,covn,':r','linewidth',2)
% hold on
for h=1:J
    semilogy(locs,covn(h,:),linestyle{h},'linewidth',2)
    hold on
end
plot(m*ones(100,1),linspace(1e-100,1),'--','color',[192,192,192]/255,'linewidth',2)
axis([-5 25 10^-2.5 1])
ylabel('C.V.')
xlabel(['Z_{',num2str(N),'}'])
hold off


% L2diff=norm(1/binwidth*Sample-exact(locs,10*myvar));
% fprintf('numsamples=%3.0f with L2diff=%8.5e \n',numSamples,L2diff)

timeElapsed=toc;
fprintf('This run took %2.2f minutes \n',timeElapsed/60)
