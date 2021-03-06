function Multiple_ISGaussian
% multiple importance sampling for the Gaussian random walk
% produces a full histogram which closely matches the exact solution
% by using multiple importance sampling targets.

% Written by Nathan Sanford to accompany Importance_Sampling_Tutorial.pdf
% The MIS sum formulas use the balance heuristic and aren't as transparent
% as those in ISGaussian.m - start there first.

% importance sampling on 
Perform_IS = 1;
% number of distributions to use
J = 6;

% details about gaussian process and numerical parameters
Ovar=1; % variance of original distribution
N=10; % number of gaussians in sum
gaussian=@(x,mu,v) 1/sqrt(2*pi*v).*exp(-(x-mu).^2/(2*v)); 
numSamples=10^5; 
numBins=101;
% locs are bin centers for histogram
locs=linspace(-25,25,numBins);
binwidth = locs(2)-locs(1);

if Perform_IS
    IS_target=linspace(-15,15,J); % biasing target(s)
    mu=IS_target/N; % mean shifts for each target
    % apportion total samples up amongst the targets
    M=floor(numSamples/J)*ones(J,1); 
    M(ceil(J/2))=M(ceil(J/2))+numSamples-sum(M);
    % pre-allocate arrays for prob and variance histograms (for each)
    PhatJ=zeros(J,length(locs)); 
    ShatJ=zeros(J,length(locs)); 
end

%% direct monte-carlo sims
binmc = zeros(1,numSamples);
% Straight monte-carlo 
for i=1:numSamples
    xmc = sqrt(Ovar)*randn(N,1);
    zz = sum(xmc);
    [~,binmc(i)] = min(abs(locs - zz));
end
% compute variance along whole distribution (there are faster ways to do
% this but this code is good for demonstration purposes)
mcP = zeros(1,numBins); mcV = mcP; 
for j = 1:numBins
    inds = find(binmc == j);
    mcP(j) = length(inds)/numSamples/binwidth;
    mcV(j) = sum(((binmc==j)-mcP(j)).^2);
end

mcV = 1/(numSamples*(numSamples-1))*mcV;
mcCV = sqrt(mcV)./mcP;

%% Using Importance sampling
if Perform_IS
    for j=1:J % this loop is for multiple importance sampling
        bin=zeros(M(j),1);
        binweight=bin;
        for k=1:M(j) % loop over samples
            % draw sample
            x=sqrt(Ovar)*randn(N,1);  % steps     
            xstar=x+mu(j); % steps with biasing
            z=sum(xstar); % random-walk with biased steps

            % binning
            finder=abs(locs-z);
            [~,bin(k)]=min(finder); 
            % likelihood ratio-each x is iid gaussian
            biased=zeros(1,J); 
            for jj=1:J
               biased(jj)=prod(gaussian(xstar,mu(jj),Ovar),1);
            end
            invlike=biased./prod(gaussian(xstar,0,Ovar)); 
            binweight(k)=1/sum(M'.*invlike);
        end

        for k=1:M(j)
            % Now do the importance sampling sum
            PhatJ(j,bin(k))=PhatJ(j,bin(k))+binweight(k); 
        end
        PhatJ(j,:)=PhatJ(j,:)./binwidth;
        for k=1:M(j)
            thisbin=bin(k);
            temp=PhatJ(j,:);
            mybinweight=M(j)*binweight(k);
            temp(thisbin)=PhatJ(j,thisbin)-mybinweight;
            ShatJ(j,:)=ShatJ(j,:)+temp.^2;
        end
    end

    Phat_IS=sum(PhatJ,1); 
    var_IS=sum(1./(M.*(M-1)).*ShatJ,1); 
    cov_IS=sqrt(var_IS)./Phat_IS; 
end

%% Plotting
if Perform_IS
    xbds=[-20 20];
    ybds=[1e-10 1];
else
    xbds=[-15 15];
    ybds=[1/(numSamples*1000) 1];
end
figure(1)
subplot(2,1,1)
semilogy(locs,gaussian(locs,0,N*Ovar),'-g',locs,mcP,'--b','linewidth',2)
legendinfo={'Exact','MC'};
hold on
if Perform_IS
    semilogy(locs,Phat_IS,':r','linewidth',2)
    for h=1:J
        semilogy(locs,PhatJ(h,:),'k')
        if h==1 || h==J
            plot(IS_target(h)*ones(100,1),linspace(1e-100,1),'--',...
                'color',0.7*ones(1,3),'linewidth',2)
        end
    end
    legendinfo{3}='MISMC estimate';
    legendinfo{4}='Individual contributions';
end
yticks(10.^(-100:2:0))
title(sprintf('10^{%d} samples of Gaussian random walk',...
    log(numSamples)/log(10)))
ylabel('P.D.F.')
axis([xbds,ybds])
hold off

subplot(2,1,2)
semilogy(1e5,1e5,'g','linewidth',2)  % for legend
hold on
semilogy(locs,mcCV,'--b','linewidth',2)
if Perform_IS
    semilogy(locs,cov_IS,':r','linewidth',2)
    plot(1e5,1e5,'k') % for legend
    for h=[1,J]
        plot(IS_target(h)*ones(100,1),linspace(1e-10,1),'--',...
            'color',0.7*ones(1,3),'linewidth',2)
    end
end
axis([xbds, min(10^-2.5,min(mcCV(mcCV>0))) 1])
% put legend in between subplots for maximal readability
myleg=legend(legendinfo,'location','north','fontsize',10.5);
myleg.Position(2) = myleg.Position(2)+0.1;  
ylabel('C.V.')
xlabel(['Z_{',num2str(N),'}'])
hold off

end
