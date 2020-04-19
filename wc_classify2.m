function outclass = wc_classify2(sample, training, group)

carry_on=true;
outidxes=1:size(sample,1);
unique_groups=unique(group);
outclass=zeros(size(outidxes));
a=0;
while carry_on
    carry_on=false;
    group_prevalence=hist(group,unique_groups);
    samples_to_classify=[];
    groups_to_classify=[];
    for g=unique(group)
        stepsize_per_group(g,:)=std(training(group==g,:),1,1);
    end
    stepsize=stepsize_per_group(group,:);
    for s=1:size(sample,1)
        distance=abs(bsxfun(@minus,training,sample(s,:)));
        distance=distance-stepsize;
        within_distance=distance<0;
        idx_close=all(within_distance,2);
        if any(idx_close)
            close_groups=group(idx_close);
            hh=hist(close_groups,unique_groups)./group_prevalence;
            [~,g]=max(hh);
            samples_to_classify=[samples_to_classify;s];
            groups_to_classify=[groups_to_classify;g];
            carry_on=true;
        end
    end
    outclass(outidxes(samples_to_classify))=groups_to_classify;
    training=[training;sample(samples_to_classify,:)];
    group=[group groups_to_classify'];
    outidxes(samples_to_classify)=[];
    sample(samples_to_classify,:)=[];
end

% 
% if nargin < 3
%     error(message('stats:classify:TooFewInputs'));
% end
% 
% % grp2idx sorts a numeric grouping var ascending, and a string grouping
% % var by order of first occurrence
% [gindex,groups,glevels] = grp2idx(group);
% nans = find(isnan(gindex));
% if ~isempty(nans)
%     training(nans,:) = [];
%     gindex(nans) = [];
% end
% ngroups = length(groups);
% gsize = hist(gindex,1:ngroups);
% nonemptygroups = find(gsize>0);
% nusedgroups = length(nonemptygroups);
% if ngroups > nusedgroups
%     warning(message('stats:classify:EmptyGroups'));
% 
% end
% 
% [n,d] = size(training);
% if size(gindex,1) ~= n
%     error(message('stats:classify:TrGrpSizeMismatch'));
% elseif isempty(sample)
%     sample = zeros(0,d,class(sample));  % accept any empty array but force correct size
% elseif size(sample,2) ~= d
%     error(message('stats:classify:SampleTrColSizeMismatch'));
% end
% m = size(sample,1);
% 
% if nargin < 4 || isempty(type)
%     type = 'linear';
% elseif ischar(type)
%     types = {'linear','quadratic','diaglinear','diagquadratic','mahalanobis'};
%     type = internal.stats.getParamVal(type,types,'TYPE');
% else
%     error(message('stats:classify:BadType'));
% end
% 
% % Default to a uniform prior
% if nargin < 5 || isempty(prior)
%     prior = ones(1, ngroups) / nusedgroups;
%     prior(gsize==0) = 0;
%     % Estimate prior from relative group sizes
% elseif ischar(prior) && strncmpi(prior,'empirical',length(prior))
%     %~isempty(strmatch(lower(prior), 'empirical'))
%     prior = gsize(:)' / sum(gsize);
%     % Explicit prior
% elseif isnumeric(prior)
%     if min(size(prior)) ~= 1 || max(size(prior)) ~= ngroups
%         error(message('stats:classify:GrpPriorSizeMismatch'));
%     elseif any(prior < 0)
%         error(message('stats:classify:BadPrior'));
%     end
%     %drop empty groups
%     prior(gsize==0)=0;
%     prior = prior(:)' / sum(prior); % force a normalized row vector
% elseif isstruct(prior)
%     [pgindex,pgroups] = grp2idx(prior.group);
%    
%     ord = NaN(1,ngroups);
%     for i = 1:ngroups
%       j = find(strcmp(groups(i), pgroups(pgindex)));
%         if ~isempty(j)
%             ord(i) = j;
%         end
%     end
%     if any(isnan(ord))
%         error(message('stats:classify:PriorBadGrpup'));
%     end
%     prior = prior.prob(ord);
%     if any(prior < 0)
%         error(message('stats:classify:PriorBadProb'));
%     end
%     prior(gsize==0)=0;
%     prior = prior(:)' / sum(prior); % force a normalized row vector
% else
%     error(message('stats:classify:BadPriorType'));
% end
% 
% % Add training data to sample for error rate estimation
% if nargout > 1
%     sample = [sample; training];
%     mm = m+n;
% else
%     mm = m;
% end
% 
% gmeans = NaN(ngroups, d);
% for k = nonemptygroups
%     gmeans(k,:) = mean(training(gindex==k,:),1);
% end
% 
% D = NaN(mm, ngroups);
% isquadratic = false;
% switch type
%     case 'linear'
%         if n <= nusedgroups
%             error(message('stats:classify:NTrainingTooSmall'));
%         end
%         % Pooled estimate of covariance.  Do not do pivoting, so that A can be
%         % computed without unpermuting.  Instead use SVD to find rank of R.
%         [Q,R] = qr(training - gmeans(gindex,:), 0);
%         R = R / sqrt(n - nusedgroups); % SigmaHat = R'*R
%         s = svd(R);
%         if any(s <= max(n,d) * eps(max(s)))
%             a=1;%error(message('stats:classify:BadLinearVar'));
%         end
%         logDetSigma = 2*sum(log(s)); % avoid over/underflow
%         
%         % MVN relative log posterior density, by group, for each sample
%         for k = nonemptygroups
%             A = bsxfun(@minus,sample, gmeans(k,:)) / R;
%             D(:,k) = log(prior(k)) - .5*(sum(A .* A, 2) + logDetSigma);
%         end
%         
%     case 'diaglinear'
%         if n <= nusedgroups
%             error(message('stats:classify:NTrainingTooSmall'));
%         end
%         % Pooled estimate of variance: SigmaHat = diag(S.^2)
%         S = std(training - gmeans(gindex,:)) * sqrt((n-1)./(n-nusedgroups));
%         
%         if any(S <= n * eps(max(S)))
%             error(message('stats:classify:BadDiagLinearVar'));
%         end
%         logDetSigma = 2*sum(log(S)); % avoid over/underflow
%         
%         if nargout >= 5
%             R = S';
%         end
%         
%         % MVN relative log posterior density, by group, for each sample
%         for k = nonemptygroups
%             A=bsxfun(@times, bsxfun(@minus,sample,gmeans(k,:)),1./S);
%             D(:,k) = log(prior(k)) - .5*(sum(A .* A, 2) + logDetSigma);
%         end
%         
%     case {'quadratic' 'mahalanobis'}
%         if any(gsize == 1)
%             error(message('stats:classify:BadTraining'));
%         end
%         isquadratic = true;
%         logDetSigma = zeros(ngroups,1,class(training));
%         if nargout >= 5
%           R = zeros(d,d,ngroups,class(training));
%         end
%         for k = nonemptygroups
%             % Stratified estimate of covariance.  Do not do pivoting, so that A
%             % can be computed without unpermuting.  Instead use SVD to find rank
%             % of R.
%             [Q,Rk] = qr(bsxfun(@minus,training(gindex==k,:),gmeans(k,:)), 0);
%             Rk = Rk / sqrt(gsize(k) - 1); % SigmaHat = R'*R
%             s = svd(Rk);
%             if any(s <= max(gsize(k),d) * eps(max(s)))
%                 error(message('stats:classify:BadQuadVar'));
%             end
%             logDetSigma(k) = 2*sum(log(s)); % avoid over/underflow
%             
%             A = bsxfun(@minus, sample, gmeans(k,:)) /Rk;
%             switch type
%                 case 'quadratic'
%                     % MVN relative log posterior density, by group, for each sample
%                     D(:,k) = log(prior(k)) - .5*(sum(A .* A, 2) + logDetSigma(k));
%                 case 'mahalanobis'
%                     % Negative squared Mahalanobis distance, by group, for each
%                     % sample.  Prior probabilities are not used
%                     D(:,k) = -sum(A .* A, 2);
%             end
%             if nargout >=5 && ~isempty(Rk) 
%                 R(:,:,k) = inv(Rk);
%             end
%         end
%         
%     case 'diagquadratic'
%         if any(gsize == 1)
%             error(message('stats:classify:BadTraining'));
%         end
%         isquadratic = true;
%         logDetSigma = zeros(ngroups,1,class(training));
%         if nargout >= 5
%             R = zeros(d,1,ngroups,class(training));
%         end
%         for k = nonemptygroups
%             % Stratified estimate of variance:  SigmaHat = diag(S.^2)
%             S = std(training(gindex==k,:));
%             if any(S <= gsize(k) * eps(max(S)))
%                 error(message('stats:classify:BadDiagQuadVar'));
%             end
%             logDetSigma(k) = 2*sum(log(S)); % avoid over/underflow
%             
%             % MVN relative log posterior density, by group, for each sample
%             A=bsxfun(@times, bsxfun(@minus,sample,gmeans(k,:)),1./S);
%             D(:,k) = log(prior(k)) - .5*(sum(A .* A, 2) + logDetSigma(k));
%             if nargout >= 5
%               R(:,:,k) = 1./S';
%             end
%         end
% end
% 
% % find nearest group to each observation in sample data
% [maxD,outclass] = max(D, [], 2);
% 
% % Compute apparent error rate: percentage of training data that
% % are misclassified, weighted by the prior probabilities for the groups.
% if nargout > 1
%     trclass = outclass(m+(1:n));
%     outclass = outclass(1:m);
%     
%     miss = trclass ~= gindex;
%     e = zeros(ngroups,1);
%     for k = nonemptygroups
%         e(k) = sum(miss(gindex==k)) / gsize(k);
%     end
%     err = prior*e;
% end
% 
% if nargout > 2
%     if strcmp(type, 'mahalanobis')
%         % Mahalanobis discrimination does not use the densities, so it's
%         % possible that the posterior probs could disagree with the
%         % classification.
%         posterior = [];
%         logp = [];
%     else
%         % Bayes' rule: first compute p{x,G_j} = p{x|G_j}Pr{G_j} ...
%         % (scaled by max(p{x,G_j}) to avoid over/underflow)
%         P = exp(bsxfun(@minus,D(1:m,:),maxD(1:m)));
%         sumP = nansum(P,2);
%         % ... then Pr{G_j|x) = p(x,G_j} / sum(p(x,G_j}) ...
%         % (numer and denom are both scaled, so it cancels out)
%         posterior = bsxfun(@times,P,1./(sumP));
%         if nargout > 3
%             % ... and unconditional p(x) = sum(p(x,G_j}).
%             % (remove the scale factor)
%             logp = log(sumP) + maxD(1:m) - .5*d*log(2*pi);
%         end
%     end
% end
% 
% %Convert outclass back to original grouping variable type
%  outclass = glevels(outclass,:);
% 
% if nargout>=5
%     pairs = combnk(nonemptygroups,2)';
%     npairs = size(pairs,2);
%     K = zeros(1,npairs,class(training));
%     L = zeros(d,npairs,class(training));
%     if ~isquadratic
%         % Methods with equal covariances across groups
%         for j=1:npairs
%             % Compute const (K) and linear (L) coefficients for
%             % discriminating between each pair of groups
%             i1 = pairs(1,j);
%             i2 = pairs(2,j);
%             mu1 = gmeans(i1,:)';
%             mu2 = gmeans(i2,:)';
%             if ~strcmp(type,'diaglinear')
%                 b = R \ ((R') \ (mu1 - mu2));
%             else
%                 b = (1./R).^2 .*(mu1 -mu2);
%             end
%             L(:,j) = b;
%             K(j) = 0.5 * (mu1 + mu2)' * b;
%         end
%     else
%         % Methods with separate covariances for each group
%         Q = zeros(d,d,npairs,class(training));
%         for j=1:npairs
%             % As above, but compute quadratic (Q) coefficients as well
%             i1 = pairs(1,j);
%             i2 = pairs(2,j);
%             mu1 = gmeans(i1,:)';
%             mu2 = gmeans(i2,:)';
%             R1i = R(:,:,i1);    % note here the R array contains inverses
%             R2i = R(:,:,i2);
%             if ~strcmp(type,'diagquadratic')
%                 Rm1 = R1i' * mu1;
%                 Rm2 = R2i' * mu2;
%             else
%                 Rm1 = R1i .* mu1;
%                 Rm2 = R2i .* mu2;
%             end
%             K(j) = 0.5 * (sum(Rm1.^2) - sum(Rm2.^2));
%             if ~strcmp(type, 'mahalanobis')
%                 K(j) = K(j) + 0.5 * (logDetSigma(i1)-logDetSigma(i2));
%             end
%             if ~strcmp(type,'diagquadratic')
%                 L(:,j) = (R1i*Rm1 - R2i*Rm2);
%                 Q(:,:,j) = -0.5 * (R1i*R1i' - R2i*R2i');
%             else
%                 L(:,j) = (R1i.*Rm1 - R2i.*Rm2);
%                 Q(:,:,j) = -0.5 * diag(R1i.^2 - R2i.^2);
%             end
%         end
%     end
%     
%     % For all except Mahalanobis, adjust for the priors
%     if ~strcmp(type, 'mahalanobis')
%         K = K - log(prior(pairs(1,:))) + log(prior(pairs(2,:)));
%     end
%     
%     % Return information as a structure
%     coeffs = struct('type',repmat({type},ngroups,ngroups));
%     for k=1:npairs
%         i = pairs(1,k);
%         j = pairs(2,k);
%         coeffs(i,j).name1 = glevels(i,:);
%         coeffs(i,j).name2 = glevels(j,:);
%         coeffs(i,j).const = -K(k);
%         coeffs(i,j).linear = L(:,k);
%         
%         coeffs(j,i).name1 = glevels(j,:);
%         coeffs(j,i).name2 = glevels(i,:);
%         coeffs(j,i).const = K(k);
%         coeffs(j,i).linear = -L(:,k);
%         
%         if isquadratic
%             coeffs(i,j).quadratic = Q(:,:,k);
%             coeffs(j,i).quadratic = -Q(:,:,k);
%         end
%     end
end

