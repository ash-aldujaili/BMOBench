function [pf,ps,fc,c]= MODIRECT(problem,l , u, numEvaluations, k)
%% FUNCTION MODIRECT
%   A multi-objective optimizer from DIRECT
% problem : function handler problem : R^n -> R^k, gives a vector 1xk
% l : lower bound of the dimension space, 1xn
% u : upper bound of the dimension space, 1xn
% numEvaluations : evaluation budget
% k : length of objective vectors
% mode : mode of the search

%-- Initialize the variables --------------------------------------%
% problem specific 
delta    = u - l;
n        = length(l);
fCount   = 0;
func     = @(x) problem( delta.*x + l);
LEVEL_TH = 5;
% tree structure
level = zeros(numEvaluations,1); % level of the node, determines its expansion and selection
c = zeros(numEvaluations,n); % centre point of the node, representative state
fc = zeros(numEvaluations, k); % f-value at that centre 1*k
dimSeq = bsxfun(@times, 1 :  n,ones(numEvaluations, n)); % the sequence of the dimension splits which the node should follow in expansion
% other variables
I = eye(n);
pf = []; ps = [];
idx = [];
% -- Initialize --
front = 1;
c(1,:) = 0.5 * ones(1,n);
fc(1,:) = func(c(1,:));
pf(1,:) = fc(1,:);
fCount = 1;
level(1) = 0;
dimSeq(1,:) = 1:n;
verbose = false;
if verbose && n ==2
    figure(2)
    rectangle('Position', [ 0 0 1 1])
    hold on
end
% -- Main Loop ---
while (fCount < numEvaluations)
	%----------------------------------------------------------------
	% 1. expand the potentially optimal nodes, one by one:
	%----------------------------------------------------------------
	for fId = 1 : length(front)
		i = front(fId);
		% sample 
		offset = (1/3)^(1 + floor(level(i)/n)); % compute the offset of the samples from the centre point
		% book-keeping variables
		firstDim = mod(level(i),n) + 1;
		numDim =  (n - firstDim + 1);
		numSamples = 2 * numDim;
		firstIdx = fCount + 1;
		lastIdx = fCount + numSamples;
		% check if you need more space for the coming samples:
		if (lastIdx > numEvaluations)
			level(firstIdx : lastIdx) = zeros(numSamples,1);
			c(firstIdx : lastIdx,:) = zeros(numSamples,n);
			fc(firstIdx : lastIdx,:) = zeros(numSamples, k);
			dimSeq(firstIdx : lastIdx,:) = bsxfun(@times, 1 :  n,ones(numSamples, n));
           % nExpnsn(firstIdx : lastIdx) = zeros(numSamples,1);
		end
		c(firstIdx : lastIdx, :) = bsxfun(@plus, c(i, :), [offset .* I(dimSeq(i,firstDim : n),:); -offset .* I(dimSeq(i,firstDim : n),:)]);
		% evaluate samples
		for j = 1 : numDim
			fc(fCount + j,:) = func(c(fCount + j,:));
			fc(fCount + numDim + j, :)= func(c(fCount + numDim + j,:));
		end
		% partition hueristic  rule : start the split along the dimension whose samples are dominating towards dominated ones
		% a.[~,idx] = sort(sum(reshape(paretofront(fc(firstIdx : lastIdx, :)), [], 2), 2), 'descend'); 
		% b. according to the distance to the parent hyperrectangle
        [~,idx] = sort(sum(reshape(pdist2(fc(firstIdx : lastIdx, :),fc(i,:)), [], 2), 2), 'descend');
        idx = dimSeq(i,idx + firstDim - 1);
		% compute the level for the child nodes and their division sequence
		for j = 1 : numDim
			offset = find(idx == dimSeq(i,j + firstDim - 1));
			level(fCount + j,:) = level(i) + offset;
			level(fCount + numDim + j, :)= level(i) + offset;
			if (mod(level(fCount+j,:),n))
			  dimSeq(fCount + j , end - numDim + 1:end) = (idx);
			  dimSeq(fCount + numDim + j, end - numDim + 1:end) = (idx);
			end
		end
		% update the level and dimseq of i
		level(i) = level(i) + numDim;
		dimSeq(i,:) = 1:n;
		% visualize
		if (verbose && (n == 2 || n ==3))
			figure(2)
            xlim([0 1]);
            ylim([0 1]);
			set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            if firstIdx==2
                firstIdx = 1;
            end
			if (n==2) % 2 dimension
				scatter(c(firstIdx : lastIdx,1),c(firstIdx : lastIdx,2),'.k')
				for j = 1 : numDim
					rectangle('Position',[c(fCount + j,1) - 3/2.*(1/3).^(1 + floor(level(fCount + j)./n) + mod(level(fCount + j),2)*(dimSeq(fCount + j,1)==1)),...
										c(fCount + j,2) - 3/2.*(1/3).^(1 + floor(level(fCount + j)./n) + mod(level(fCount + j),2)*(dimSeq(fCount + j,1)==2)),... 
										3.*(1/3).^(1 + floor(level(fCount + j)./n) + mod(level(fCount + j),2)*(dimSeq(fCount + j,1)==1)),...
										3.*(1/3).^(1 + floor(level(fCount + j)./n) + mod(level(fCount + j),2)*(dimSeq(fCount + j,1)==2))]);
					rectangle('Position',[c(fCount + numDim + j,1) - 3/2.*(1/3).^(1 + floor(level(fCount + numDim + j)./n) + mod(level(fCount + numDim + j),2)*(dimSeq(fCount + numDim+ j,1)==1)),...
										c(fCount + numDim + j,2) - 3/2.*(1/3).^(1 + floor(level(fCount + numDim + j)./n) + mod(level(fCount + numDim + j),2)*(dimSeq(fCount + numDim+ j,1)==2)),... 
										3.*(1/3).^(1 + floor(level(fCount + numDim + j)./n) + mod(level(fCount + numDim + j),2)*(dimSeq(fCount + numDim+ j,1)==1)),...
										3.*(1/3).^(1 + floor(level(fCount + numDim + j)./n) + mod(level(fCount + numDim + j),2)*(dimSeq(fCount + numDim+ j,1)==2))]);

				end
			else % 3 dimension
				scatter3(c(firstIdx : lastIdx,1),c(firstIdx : lastIdx,2),c(firstIdx : lastIdx,3),'.k')
			end
			hold on
			
			%keyboard
			%pause(0.1);
       end
       	
		% update fCount
		fCount = fCount + numSamples;

		if (fCount >= numEvaluations)
			break;
		end
	end
	%----------------------------------------------------------------
	% 2. identify the potentially optimal nodes:
	%----------------------------------------------------------------
% 	% add attributes to the fc (pseudo-function space) and look for the non-dominated samples according to that space:
%     if (strcmp(mode, 'EE')) % explore - exploit
% 	  	%front = find(paretofront([ prod(bsxfun(@minus,fc(1:fCount,:),min(fc(1:fCount,:),[],1)-1),2) level(1:fCount)]));
%   %     pArray = [fc(1:fCount,:), nExpnsn(1:fCount), level(1:fCount) prod(bsxfun(@minus,fc(1:fCount,:),min(fc(1:fCount,:),[],1)-1),2)];
%          pArray = [level(1:fCount) ndsort(fc(1:fCount,:))];
%       % pArray(front,:) = inf;
%      % keyboard 
%        front = find(paretofront(pArray));
%     elseif (strcmp(mode, 'EED')) % explore - exploit - discover
% 	  	front = find(paretofront([ fc(1:fCount,:) prod(bsxfun(@minus,fc(1:fCount,:),min(fc(1:fCount,:),[],1)-1),2) level(1:fCount) -sum((pdist2(fc(1:fCount,:),pf)),2)]));
%     elseif (strcmp(mode,'ED'))
%         front = find(paretofront([fc(1:fCount,:) -sum((pdist2(fc(1:fCount,:),pf)),2)]));
%     end
    % select those above a size threshold

    %front = find(paretofront([fc(1:fCount,:) level(1:fCount)]));
    levelTh = min(level(1:fCount)) + LEVEL_TH;
    idx = level(1:fCount) <= levelTh;
    %keyboard
    frontIdx = (paretofront([fc(idx,:) level(idx)]));
    front = zeros(fCount,1);
    front(idx) = frontIdx;
    front = find(front);

    %front = find(paretofront([fc(1:fCount,:) level(1:fCount)]));
    % update the expansion counter for each
	%----------------------------------------------------------------
	% 3. Get the approximation set
	%frontSet = paretofront(fc);
	%pf = fc(frontSet,:);
	%----------------------------------------------------------------
            %fc(firstIdx:lastIdx,:)
          % keyboard
end
% -- End of Main Loop --
% Get the approximation set
c = bsxfun(@plus, bsxfun(@times,delta, c),l); % normalize
front = paretofront(fc);
pf = fc(front,:);
ps = c(front,:);


end
