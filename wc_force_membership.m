function class_out = wc_force_membership(f_in, class_in, f_out, handles)
% class = function force_membership(f_in, class_in, f_out, handles)
% Given classified points, try to classify new points via template matching
%
% f_in:          features of classified points  (# input spikes x n_features)
% class_in:      classification of those points
% f_out:         features of points to be classified (nspk x n_features)
% handles        environment variables, of which the following are
%                required: 
%                    o handles.par.template_sdnum - max radius of cluster,
%                                                   in std devs.
%                    o handles.par.template_k     - # of nearest neighbors
%                    o handles.par.template_k_min - min # of nn for vote
%                    o handles.par.template_type  - nn, center, ml, mahal

nspk = size(f_out,1);
class_out = zeros(1,size(f_out,1));
switch handles.par.template_type
    case 'nn'
        sdnum = handles.par.template_sdnum;
        k     = handles.par.template_k;
        k_min = handles.par.template_k_min;
        sd    = sqrt(sum(var(f_in,1)))*ones(1,size(f_in,1));
        for i=1:nspk,
            nn = nearest_neighbor(f_out(i,:),f_in,sdnum*sd,Inf*ones(size(f_in)),Inf,k);
            if( nn )
                winner = vote(class_in(nn),k_min);
                class_out(i) = winner;
            else
                class_out(i) = 0;
            end
        end
      
    case 'center'
        [centers, sd, pd] = build_templates(class_in,f_in); % we are going to ignore pd
        sdnum = handles.par.template_sdnum;
        for i=1:nspk,
            class_out(i) = nearest_neighbor(f_out(i,:),centers,sdnum*sd);        
        end
        
    case 'ml'
        [mu sigma] = fit_gaussian(f_in,class_in);
        for i=1:nspk,
            class_out(i) = ML_gaussian(f_out(i,:),mu,sigma);
        end
    case 'mahal'
        [mu sigma] = fit_gaussian(f_in,class_in);
        for i=1:nspk,
            class_out(i) = nearest_mahal(f_out(i,:),mu,sigma);
        end
        
    otherwise
        sprintf('force_membership(): <%s> is not a known template type.\n',handles.par.template_type);
        
end

function [mu, sigma] = fit_gaussian(x,class)

N = max(class);
mu = zeros(N,size(x,2));
sigma = zeros(size(x,2),size(x,2),N);

for i=1:N,
    mu(i,:) = mean(x(class==i,:));
    sigma(:,:,i) = cov(x(class==i,:));
end

function index = ML_gaussian(x,mu,sigma)
% function index = ML_gaussian(x,mu,sigma)
% x is a vector drawn from some multivariate gaussian
% mu(i,:) is the mean of the ith Gaussian
% sigma(:,:,i) is the covariance of the ith Gaussian
% 
% Returns the index of the Gaussian with the highest value of p(x).

N = size(mu,1);  % number of Gaussians

if( N == 0 )
    index = 0;
else
    for i=1:N,
        % leave out factor of 1/(2*pi)^(N/2) since it doesn't affect argmax
        p(i) = 1/sqrt(det(sigma(:,:,i)))*exp(-0.5*(x-mu(i,:))*inv(sigma(:,:,i))*(x-mu(i,:))');
    end
    [m index] = max(p);
end

function index = nearest_mahal(x,mu,sigma)
% function index = nearest_mahal(x,mu,sigma)
% x is a vector
% mu(i,:) is the mean of the ith Gaussian
% sigma(:,:,i) is the covariance of the ith Gaussian
% 
% Returns the index of the Gaussian closest (by the Mahalanobis distance)
% to x.

N = size(mu,1);  % number of Gaussians
d = [];
if( N == 0 )
    index = 0;
else
    for i=1:N,
        d(i) = (x-mu(i,:))*inv(sigma(:,:,i))*(x-mu(i,:))';
    end
    [m index] = min(d);
end

function index = nearest_neighbor(x,vectors,maxdist,varargin)
% function index = nearest_neigbor(x,vectors,maxdist,pointdist*,pointlimit*,k*)
% x is a row vector
% pointdist (optional) - vector of standard deviations
% pointlimit (optional) - upper bound on number of points outside pointdist
% k (optional) - number of points used for nearest neighbor

% Find the distance to all neighbors. Consider only those neighbors where
% the point falls in the radius of possibility for that point. Find the
% nearest possible neighbor.
% Return 0 if there is no possible nearest neighbor.

distances = sqrt(sum((ones(size(vectors,1),1)*x - vectors).^2,2)');
conforming = find(distances < maxdist);
if( length(varargin) > 0 )
    pointdist = varargin{1};
    if( length(varargin) > 1 )
        pointlimit = varargin{2};
    else
        pointlimit = Inf;
    end
    pointwise_conforming = [];
    for i=1:size(vectors,1),
        if( sum( abs(x-vectors(i,:)) > pointdist(i,:) ) < pointlimit )  % number of deviations from pointdist allowed.
            pointwise_conforming = [pointwise_conforming i];
        end
    end
    conforming = intersect(conforming, pointwise_conforming);
end
if( length( conforming ) == 0 )
    index = 0;
else
    if( length(varargin) > 2 )
        k = varargin{3};
        [y i] = sort(distances(conforming)); % k-nearest neighbors
        i = i(1:min(length(i),k));
    else
        [y i] = min(distances(conforming));   
    end
    index = conforming(i);
end


