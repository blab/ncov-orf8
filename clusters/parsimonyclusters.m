function [] = parsimonyclusters(treePath,metaPath,koCol,koValue,outFile)
    tree = phytreeread(treePath);
    % get all WA sequences
    f = fopen(metaPath);
    col = strsplit(fgets(f), '\t');
    strain_id = find(ismember(col,'strain'));
    ko_id = find(contains(col,koCol));
    c=1;
    id = cell(0,0); 
    strain = cell(0,0); 
    ko = cell(0,0); 
    while ~feof(f)
        col = strsplit(fgets(f), '\t','CollapseDelimiters', false);
        id{c,1} = col{1};
        strain{c,1} = col{strain_id};
        if koCol == 'ORF9b_koType'
            ko{c,1} = erase(col{ko_id},newline);
        else
            ko{c,1} = col{ko_id};
        end
        c=c+1;
    end
    fclose(f);
    
    % get all leafnames
    leafs = get(tree, 'leafnames');
    %%
    is_ko = find(ismember(ko,koValue));
    
    is_ko_leaf = false(length(leafs), 1);
    for i = 1 : length(leafs)
        if sum(ismember(id(is_ko),leafs{i}))==1
            is_ko_leaf(i) = true;
        end
    end
    
    nodenames = get(tree, 'nodenames');
    
    % get the connectivity matrix and the distances between nodes
    mat = getmatrix(tree);
    %dist = get(tree, 'Distances');
    %% initialize the location vector
    
    location = cell(length(nodenames),1);
    visited = false(length(nodenames),1);
    
    for j = 1 : length(is_ko_leaf)
        visited(j) = true;
        if is_ko_leaf(j)
            location{j} = "ko";
        else
            location{j} = "noko";
        end
    end
    %%
    % upwards path of the parsimony calling
    not_visited = find(~visited);
    while ~isempty(not_visited)
        disp(length(not_visited))
        not_visited = find(~visited);
        for j = 1 : length(not_visited)
            children = find(mat(not_visited(j),:));
            if ~isempty(location{children(1)}) && ~isempty(location{children(2)})
                    int = intersect(location{children(1)}, location{children(2)});
                    if isempty(int)
                        location{not_visited(j)} = [location{children(1)}, location{children(2)}];
                    else
                        location{not_visited(j)} = int;
                    end
                visited(not_visited(j)) = true;
            end
        end
    end
    %%
    %% downwards calling
    
    visited(length(is_ko_leaf)+1:end) = false;
    visited(end) = true;
    
    not_visited = find(~visited);
    while ~isempty(not_visited)
        not_visited = find(~visited);
        for j = length(not_visited) : -1 : 1
            parent = find(mat(:,not_visited(j)));
    
            uni_loc = unique(location{not_visited(j)});
    
            if sum(ismember(not_visited, parent))==0
                if length(uni_loc) < length(location{not_visited(j)})
                    combined_locs = [location{not_visited(j)} location{parent}];
                    uni_comb = unique(combined_locs);
                    freqs = zeros(length(uni_comb),1);
                    for k = 1 : length(uni_comb)
                        freqs(k) = sum(combined_locs==uni_comb{k});
                    end
                    ind_max = find(freqs==max(freqs));
                    location{not_visited(j)} = uni_comb(ind_max);
                else
                    int = intersect(location{not_visited(j)}, location{parent});
                    if ~isempty(int)
                        location{not_visited(j)} = int;
                    end
                end
                if length(location{not_visited(j)})>1
                    location{not_visited(j)};
                end
                visited(not_visited(j)) = true;
            end
        end
    end
    %%
    % get all nodes that are only orf8KOs
    onlyKO = false(length(visited),1);
    uncertainCOunt = 0;
    for j = 1 : length(onlyKO)
        if length(location{j})==1 && location{j}=="ko"
            onlyKO(j) = true;
        elseif length(location{j})==2
            onlyKO(j) = true;
            uncertainCOunt = uncertainCOunt+1;
        end
    end
    
    %%
    isInKO = find(onlyKO(1:length(is_ko_leaf)));
    
    ko_names = str2double(nodenames(isInKO));
    
    clustered_node = false(size(nodenames));
    clustered_node(isInKO) = true;
    
    
    % for each isInWA leaf, get all parent nodes in WA
    KOParent = cell(length(isInKO),1);
    for a = 1 : length(isInKO)
        parent = find(mat(:,isInKO(a)));
        while onlyKO(parent)
            clustered_node(parent) = true;
            KOParent{a} = [KOParent{a} parent];
            parent = find(mat(:,parent));
        end
    end
    
    in_cluster_mat = zeros(length(isInKO), length(isInKO));
    
    new_mat = zeros(size(in_cluster_mat));
    % get for each pair of leaves if they are in the same cluster
    for a = 1 : length(isInKO)
        %disp(a)
        if ~isempty(KOParent{a})
            for b = a+1 : length(isInKO)
                if ~isempty(KOParent{b})
                    int = intersect(KOParent{a}, KOParent{b});
                    if ~isempty(int)
                        new_mat(a,b) = 1;
                    end
                end
            end
        end
    end
    %%
    ko_leafs = leafs(is_ko_leaf);
    
    connected = new_mat;
    
    stop=false;
    while ~stop
        [a,b] = find(connected==1,1,'first');
        if isempty(a)
            stop = true;
        else
            % combine a and b
            ko_leafs{a} = [ko_leafs{a} ',' ko_leafs{b}];
            % combine the rows of a and b
            for i = 1 : length(ko_leafs)
                connected(a,i) = max([connected(a,i),connected(b,i)]);
            end
            % remove b
            connected(b,:) = 0;
            connected(:,b) = 0;
            ko_leafs{b} = 'NA';
        end
    end
    
    % get the cluster
    cl_ind = find(~ismember(ko_leafs, 'NA'));
    ko_clusters = ko_leafs(cl_ind);
    
    
    f = fopen(outFile, 'w');
    fprintf(f,'strain\tcluster\n');
    for a = 1 : length(ko_clusters)
        seqs = strsplit(ko_clusters{a}, ',');
        for b = 1 : length(seqs)
            if length(seqs{b})>2
                fprintf(f, '%s\t%d\n', seqs{b}, a);
            end
        end
    end
    fclose(f);
end