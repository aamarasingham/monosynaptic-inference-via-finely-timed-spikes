function [CIs] = cisynch(r,t,delta,alpha,maxiter)
% Calculate assumption-less confidence intervals ("worst case").
% This function assumes your spike times are rounded to milliseconds.
% Also, please note that this function defines "synchrony" as spikes
% times that are exactly equal. For lagged sychrony you must send 
% the function a lagged version of t (it is pretty fast).
% Lastly, please beware, this function might not be numerically stable for
% very large spike trains. However, the computation time technically does not
% depend on the duration of the the spike trains, but rather on #
% the co-variance of spiking density.
% For further questions contact:
% Zach Saccomano, zsaccomano4@gmail.com
% ------ INPUTS --------
%   'r'    -  round reference spike times (putative pre-synaptic spike times)
%   't'    -  round target spike times (putative post-synaptic spike times)
%   'delta'   -  background timescale (in milliseconds)
%   of milliseconds).
%   'alpha'   -  Confidence level (i.e. 0.05 gives you two-sided 95%
%   confidence intervals.
%   'maxiter' - maximal number of iterations for optimization process. In
%   reality, I find that the exact solution converges iterations<20... but I
%   still put this in for the sake of formality.

% ------ OUTPUTS --------
%   'CIs'   -  1 minus alpha level confidence intervals (in units of number of synchrony).


S = numel(intersect(r,t)); % total observed synchrony
bEdges = 0:(delta):delta*floor(max(t)./delta);
refCounts = histc(r,bEdges); % reference counts
w = floor(t./delta) + 1; % window
Nr = refCounts(w); % reference counts for target spike i
Nr = Nr(Nr~=0); % they don't contribute
p = sort(Nr,'descend')./delta; % bernoulli weights
q = flipud(p); % for second tail
uniqP = unique(p)';
M = [uniqP;1-uniqP;zeros(numel(p),numel(uniqP))]; % pre-allocate FFT matrix.
df = fft(M); % only need to take FFT once for the unique p's

% Compute CIs
for b = 1:2 % for lower and upper bound
    thta = 1; val = 0; sePnt = 1; D = 2; prevThta = thta; % set optimization points
    for k = 1:maxiter % don't optimize past maxiter partitions
        thta = round(length(p).*(sePnt./D)); % test value of theta
        if thta==0
            thetaFin(b) = 0;
           break;
        end
        if b == 1
            counts = histc(q(thta:end),uniqP)'; % number of redundant r.v.'s to convolve
        else
            counts = histc(p(thta:end),uniqP)'; % same for upper bound
        end
        dfP = prod(df.^repmat(counts,length(df),1),2); % scale the unique r.v.'s in the frequency domain
        dist1 = ifft(dfP); % get the distribution
        dist1 = flipud(dist1); % flip it (convention)
        x = [1:length(dist1)]; % define synchrony axis
        if b==1
            val = sum(dist1(x' >= round(S))); % p-value for lower bound
        else
            val = sum(dist1(x' <= round(S))); % same for upper bound
            
        end
        if abs(thta - prevThta) == 1 % ...if the previous iteration only changed by one theta
            if (val - alpha/2) > 0 && (prevVal - alpha/2) < 0 % ...and the value of the integral crossed alpha/2
                thetaFin(b) = prevThta; % ...the optimization is finished.
                break;
            elseif prevVal - alpha/2 > 0 && (val - alpha/2) < 0 % but it's also true in the other direction.
                thetaFin(b) =  thta;
                break;
            end
        end
        D = D.*2; % divide the search space in half every iteration
        if b == 1
            if val   <= alpha/2 % decide which way to increment down the tree structure of a geometric series
                sePnt = sePnt*2 + 1;
            elseif val >= alpha/2
                sePnt = sePnt*2 - 1;
            end
            
        else
            if val <= alpha/2
                sePnt = sePnt*2 - 1;
            elseif val >= alpha/2
                sePnt = sePnt*2 + 1;
            end
        end
        prevThta = thta; % save thta from last iteration so we know when we are done.
        prevVal = val; % and the value of the integral.
    end
end
CIs = thetaFin; % injected synchrony CI's.
end

