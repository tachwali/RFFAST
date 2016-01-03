%*********************************************************************%
%% R-FFAST: Fast Fourier Aliasing Sparse Transform Demo
% File: RFFAST.m
% Author: Yahia Tachwali, Ray Duran
% Description: This algorithm follows the techniques as outlined
% in the  paper: " A robust sub-linear time R-FFAST algoriothm for 
% computing a sparse DFT"
% by Sameer Pawart and Kannan Ramchandran, Jan 2015
% Dept of Electrical Engineering U.C. Berkeley
% 
%%**********************************************************************%
% Revision:
% Rev 0.1 rbd 12-12-15          Initial draft (assumed downsampling for three stages 3, 4 , 5)
% Rev 0.2                       Corrected versioning and added input check
%**********************************************************************%

function Xout = RFFAST(x)
    % input variable must be a vector    
    [m,n]=size(x);
    N = m;
    if( ~isvector(x) )
        error('input variable is not a vector');    
    end    
    if( m < n )
        x = x';
        N = n;
    end    
    if( rem(N,60) > 0 )
        error('Current implementation supports input variable with length that is a multiple of 60');
    end
        
    % TODO: temporary resizing of input signal to fit the assumption    
    Xout = zeros(N,1);    
    ITERATIONS = 10;
    T = 0.0001; %threshold
    NUMOFSTAGES = 3;
    NUMOFCLUSTERS = 3;
    CLUSTERSIZE = 5;
    NUMOFBRANCHES = NUMOFCLUSTERS * CLUSTERSIZE;
        
    downSamplingFactors = findDownSamplingFactors(N,NUMOFSTAGES);
    NUMOFSTAGES = length(downSamplingFactors);% update NUMOFSTAGES in case the function above return more stages based on N

    % Front-end
    for s = 1 : NUMOFSTAGES        
        %setup        
        [xs, startIndeces] = samplingStage(x, NUMOFCLUSTERS, CLUSTERSIZE, downSamplingFactors(s));
        XS_TEMP = zeros( size(xs) );
        
        % obtain FFTs
        for b = 1:NUMOFBRANCHES
            XS_TEMP(b,:) = fft(xs(b,:));
        end
        
        NUMOFBINS(s) = size(XS_TEMP, 2);
        % prepare results in terms of stages and bins
        for bin = 1 : NUMOFBINS(s)
            XS{s,bin} = XS_TEMP(:,bin);
        end
        
        AllStartIndeces{s} = startIndeces;
    end
    
    % Back-end (peeling engine)
    Xfound = zeros(N,1);
    multiToneHit = 0;%for debugging only
    
    for i = 1 : ITERATIONS
        for s = 1 : NUMOFSTAGES            
            for b = 1 : NUMOFBINS(s)
                y = XS{s,b};
                if( norm(y)^2 < T )
                    Xfound(b:downSamplingFactors(s):end)=0;
                else
                    [isSingleton, vp, p] = SingletonEstimator(y, NUMOFCLUSTERS, CLUSTERSIZE, N, AllStartIndeces{s});
                    if(isSingleton)
                        a = steering_vector(p,AllStartIndeces{s}, CLUSTERSIZE, N);
                        XS = Peel(XS, vp, p, a);
                        Xout(p) = vp;
                    else
                        multiToneHit = multiToneHit+1; %for debugging only
                    end
                end    
            end        
        end
    end        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% samplingStage
% returns xs[ numOfBranches x numOfbins ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xs, startIndeces]=samplingStage(  x, ...
                                            numOfClusters, ...
                                            clusterSize, ...
                                            downSamplingFactor)
    N=length(x);
    numOfBranches = numOfClusters * clusterSize;    
    
    xs = zeros( numOfBranches, N/downSamplingFactor );
        
    startIndeces=randi(N,1,numOfClusters);
    
    k = 1;
    for i = 1 : numOfClusters
        start_index = startIndeces(i);
        for j = 0 : (clusterSize-1)
            % even if start_index is larger than x size, circshift does the
            % right thing
            x_circularShifted = circshift(x,start_index); 
            xs(k,:)=x_circularShifted(1:downSamplingFactor:end);
            
            %update indeces
            k = k + 1;
            start_index = start_index + 2^(i-1);                        
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% findDownSamplingFactors returns a list of down-sampling factors to be
% used for each stage
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function factors = findDownSamplingFactors(N,NumberforStages)
%TODO : find factors based on N and number of stages(if necessary)
    factors = [3,4,5];
end

function [isSingleton, vp, p] = SingletonEstimator(Y, C, N, Nsamp, start_indeces)
    wprev = 0;
    isSingleton = false;    
    D = N*C;
%TODO : calculating Beta from the single frequency estimation paper
    Beta = ones(N-1,1);
    Beta = Beta/length(Beta);
    gamma = 0.1;
        
    for i = 0 : C-1
    % We will assume that we have generated phase vectors??
        wnew = 0;
        for t = (N*i)+1 : (N*(i+1)-2)+1 % Because of silly matlab :(
           	wnew = wnew + Beta(t-N*i)*angle(Y(t+1)/Y(t)) ;
        end
        
        if(wnew<0)
            wnew = 2*pi + wnew;
        end
        wnew = (1/(2^i))*wnew;
        delta_1 = ceil(wprev/((2*pi)/(2^i)))*((2*pi)/(2^i)) + wnew - wprev;
        delta_2 = floor(wprev/((2*pi)/(2^i)))*((2*pi)/(2^i)) + wnew - wprev;        
        
        if( abs(delta_1) < abs(delta_2)) 		
            wprev = delta_1 + wprev;
        else
            wprev = delta_2 + wprev;
        end
    end
	
	% Set support estimate for q
    % TODO: I got the result correct with 4pi not 2pi
	q = round(wprev*Nsamp/(2*pi)); % add abs value ??
    
    q = q+ 1; % ??? for zero
	
	%Set the energy threshold
	T = ( 1 + gamma)*D;
	
    % Calculate a(q) vectors....for specific bin
    [aq] = steering_vector(q,start_indeces, N, Nsamp);    
    % Calculate amplitude of basis functions
    vq = (conj(aq)*Y)/D;  %*' get rid of grey!
    % Take the norm
    %%X = Y - vq*aq';
    %%if ( (norm(X))^2 < T)
    X=abs(Y);
    if ( (norm(X-mean(X)))^2 < T)
    	isSingleton = true;
    	p = q;
    	vp = vq;
    else
        p = 0;
        vp = 0;
    end
end

function XS_peeled = Peel(XS,vp,p,a)
    [numOfStage, numOfBins] = size(XS);        
    for s = 1 : numOfStage                
        for b = 1 : numOfBins
            if(~isempty(XS{s,b}))
                XS{s,b} = XS{s,b} - vp*a';
            end
        end
    end    
    XS_peeled = XS;
end


function [aq] = steering_vector(q,start_indeces, clusterSize, NSamp)
    
    v = zeros(1,length(start_indeces)*clusterSize);
    k=0;
    for i = 1 : length(start_indeces)           
        for j = 0 : (clusterSize-1)            
            k = k + 1;       
            v(k) = start_indeces(i) + j*2^(i-1);               
        end
    end

    aq = exp(2*pi*sqrt(-1)*q.*v/NSamp);       
end