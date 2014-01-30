

function [Results, Num_C, Weights, Clusters, Centers] = RSK_Means(data, k, alpha, lass)




[d, n] = size(data);

w(1:d) = 1/sqrt(d); % Feature weights
C = zeros(1,n);  % Cluster assignments

stream = RandStream('mt19937ar','Seed',0);
rperm = randperm(stream, n, k);
% data = [data1 data2 data3];
cen = data(:, rperm);  % Cluster centers
num_trimmed = floor(alpha * size(data,2));


old_w = zeros(1,d);
it_cntr = 0;
% while(unweighted_dist(old_w, w) > 10^-3)
while(abs(sum(w - old_w))/sum(abs(old_w)) > 10^-3)
    it_cntr = it_cntr + 1
%     unweighted_dist(old_w, w)
    abs(sum(w - old_w))/sum(abs(old_w))
%     abs(sum(abs(w)) - sum(abs(old_w)))/sum(abs(old_w))
    old_w = w;

    %------------------------------- Step(a/b) -------------------------------
    % Perform Trimmed K-Means on the weighted dataset

    old_C = ones(size(C));
    while(unweighted_dist(old_C, C) > 0)
        old_C = C;
        % Assign each point to the cluster center with minimum weighted distance
        wdists_2cntrs = weighted_dist(w, data', cen');
        [~,C_idxs] = min(wdists_2cntrs,[],2);
        C = C_idxs';

        % Find number of samples assigned to each cluster
        number_samps = zeros(k,1);
        for i=1:k
            number_samps(i,1) = size(find(C==i),2);
        end

        % Sort the points by weighted distance to assigned cluster center
        sorted_idxs = zeros(k,max(number_samps));
        sorted_dists = zeros(k,max(number_samps));
        for i=1:k
            samp_idxs = find(C==i);
            weighted_distances = zeros(size(samp_idxs));
            for j=1:size(samp_idxs,2)
                weighted_distances(j) = weighted_dist(w, cen(:,i)', data(:,samp_idxs(j))');
            end
            mat = [samp_idxs; weighted_distances];
            [sdists, sidxs] = sort(mat(2,:),2,'descend');
            sorted_dists(i,1:number_samps(i)) = sdists;
            sorted_idxs(i,1:number_samps(i)) = samp_idxs(sidxs);
        end

        % Find the points to be trimmed
        Ow = zeros(1,num_trimmed);
        cluster_trim_idx = ones(1,k);
        for i=1:num_trimmed
            vals = zeros(1,k);
            for j=1:k
                vals(1,j) = sorted_dists(j,cluster_trim_idx(1,j));
            end
            [~, max_idx] = max(vals);
            Ow(1,i) = sorted_idxs(max_idx, cluster_trim_idx(1,max_idx));
            cluster_trim_idx(1,max_idx) = cluster_trim_idx(1,max_idx) + 1;
        end

        % Find cluster centers from non-trimmed points
        for i=1:k
            cen(:,i) = mean(data(:,sorted_idxs(i,cluster_trim_idx(1,i):number_samps(i,1))),2);
        end
    end

    %-------------------------------------------------------------------------


    %-------------------------------- Step(c) --------------------------------
    % Find points with largest unweighted distance from their cluster centers

        % Sort the points by unweighted distance to assigned cluster center
        uw_sorted_idxs = zeros(k,max(number_samps));
        uw_sorted_dists = zeros(k,max(number_samps));
        for i=1:k
            samp_idxs = find(C==i);
            unweighted_distances = zeros(size(samp_idxs));
            for j=1:size(samp_idxs,2)
                unweighted_distances(j) = unweighted_dist(cen(:,i)', data(:,samp_idxs(j))');
            end
            mat = [samp_idxs; unweighted_distances];
            [sdists, sidxs] = sort(mat(2,:),2,'descend');
            uw_sorted_dists(i,1:number_samps(i)) = sdists;
            uw_sorted_idxs(i,1:number_samps(i)) = samp_idxs(sidxs);
        end

        % Find the points to be trimmed
        Oe = zeros(1,num_trimmed);
        cluster_trim_idx = ones(1,k);
        for i=1:num_trimmed
            vals = zeros(1,k);
            for j=1:k
                vals(1,j) = uw_sorted_dists(j,cluster_trim_idx(1,j));
            end
            [max_val, max_idx] = max(vals);
            Oe(1,i) = uw_sorted_idxs(max_idx, cluster_trim_idx(1,max_idx));
            cluster_trim_idx(1,max_idx) = cluster_trim_idx(1,max_idx) + 1;
        end

    %-------------------------------------------------------------------------


    %-------------------------------- Step(d) --------------------------------
    % Combine sets Ow and Oe to get the trimmed and untrimmed points

    O = union(Ow, Oe);
    num_untrimmed = zeros(1,k);
    untrimmed_idxs = zeros(size(sorted_idxs));
    for i=1:k
        dif = setdiff(sorted_idxs(i,1:number_samps(i,1)), O);
        num_untrimmed(1,i) = size(dif,2);
        untrimmed_idxs(i,1:num_untrimmed(1,i)) = dif;
    end

    %-------------------------------------------------------------------------



    %-------------------------------- Step(e) --------------------------------
    % Solve for the weights by optimizing the BCSS

    BCSS = zeros(1,k);
    tot_untrimmed = sum(num_untrimmed);
    for i=1:d
%         weight_solve_cntr = i
        all_points = zeros(tot_untrimmed, d);
        ap_cntr = 1;
        for j=1:k
            for l=1:num_untrimmed(1,j)
                all_points(ap_cntr,:) = data(i,untrimmed_idxs(j,l))';
                ap_cntr = ap_cntr + 1;
            end
        end
        AP_Dists = unweighted_dist(all_points, all_points);
        term_one = sum(AP_Dists(:))/tot_untrimmed;
        
        term_two = 0;
        for j=1:k
            cluster_points = zeros(num_untrimmed(1,j),d);
            cp_cntr = 1;
            for l=1:num_untrimmed(1,j)
                cluster_points(cp_cntr,:) = data(i,untrimmed_idxs(j,l))';
                cp_cntr = cp_cntr + 1;
            end
            CP_Dists = unweighted_dist(cluster_points,cluster_points);
            term_two = term_two + sum(CP_Dists(:))/num_untrimmed(1,j);
        end
        
        BCSS(1,i) = term_one - term_two;
    end


    wofzero = BCSS/sqrt(BCSS*BCSS');
    if(sum(abs(wofzero)) <= lass)
        w= wofzero;
        sum(abs(wofzero))
    else
        del = get_delta(BCSS, d, lass);
        w = abs(my_plus(BCSS-del*ones(1,d))/sqrt(my_plus(BCSS-del*ones(1,d))*my_plus(BCSS-del*ones(1,d))'));
    end

    %-------------------------------------------------------------------------

end

Num_C = number_samps;
Weights = w;
Results = zeros(d,max(number_samps),k);
Clusters = C;
Centers = cen;
for i=1:k
    Results(:,1:number_samps(i,1),i) = data(:,find(C==i));
end



end