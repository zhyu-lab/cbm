function visualize_results(cellAssignFile,genotypeFile,outputFile)

fid = fopen(cellAssignFile,'r');
line = fgetl(fid);
labels_p = str2double(regexp(line,'\s','split'))+1;
fclose(fid);

clone_genotypes = load(genotypeFile);
[num_clone,num_muta] = size(clone_genotypes);
predicted = clone_genotypes(labels_p,:);

% construct MST
dists = ones(num_clone,num_clone)*10000;
for i = 1:num_clone-1
    for j = i+1:num_clone
        tmp = clone_genotypes(i,:)-clone_genotypes(j,:);
        dists(i,j) = sqrt(tmp*tmp');
        dists(j,i) = dists(i,j);
    end
end
r_indx = 1;
for i = 2:num_clone
    tmp = sum(clone_genotypes(i,:));
    if tmp < sum(clone_genotypes(r_indx,:))
        r_indx = i;
    end
end

count = 1;
added = r_indx;
all = 1:num_clone;
trees = zeros(1, num_clone);
while count < num_clone
    tv = ismember(all,added);
    no_added = all(~tv);
    tmp = dists(added,no_added);
    [~, I] = min(tmp(:));
    i = rem(I-1,length(added))+1;
    j = floor((I-1)/length(added))+1;
    added = [added no_added(j)];
    trees(no_added(j)) = added(i);
    count = count+1;
end

[coef,score,latent] = pca(clone_genotypes');
tmp = cumsum(latent);
props = tmp/tmp(end);
indxs = find(props > 0.9);
I = indxs(1);
score = score(:,1:I);
I = kmeans(score,num_clone,'EmptyAction','singleton','Replicates',200);
sorted_muta_indxs = [];
muta_indxs = 1:num_muta;
for i = 1:num_clone
    tv = I == i;
    sorted_muta_indxs = [sorted_muta_indxs muta_indxs(tv)];
end

sorted_cell_indxs = [];
for i = added(end:-1:1)
    tv = labels_p == i;
    sorted_cell_indxs = [sorted_cell_indxs find(tv)];
end

fontsize = 18;
predicted = predicted(sorted_cell_indxs,sorted_muta_indxs);
h1 = HeatMap(predicted, 'Colormap', colormap([1 1 1; 0 0 0]), 'Symmetric', 'false');
addXLabel(h1, 'Mutation', 'fontsize', fontsize);
addYLabel(h1, 'Cell', 'fontsize', fontsize);

fh = figure(1);
plot(h1,fh);

print(fh,outputFile,'-dpng','-r300')

close all;

end