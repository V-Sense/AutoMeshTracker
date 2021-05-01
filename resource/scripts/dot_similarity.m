
function keyframe_indices = dot_similarity(descriptor_path,feas_path,max_region_length)

    % grab the output path
    [filepath] = fileparts(char(descriptor_path));

    descriptors = load(descriptor_path);

    for row_idx = 1:size(descriptors,1)
        for row_idx_2 = 1:size(descriptors,1)
        
            %Generate similarity using simple dot product of features
            similarity(row_idx,row_idx_2)=dot(descriptors(row_idx,:),descriptors(row_idx_2,:));
            
        end
    end
    
    % normalize and plot
    sim_norm=normalize(similarity,'range');
    heatmap(sim_norm);
    
    % get average sim score per frame
    averaged_signal = mean(sim_norm,2);
    
    % load feas score
    feas = load(feas_path);

    % filter data with moving average
    % define coeffs
    a=1;
    % b=repmat(1/4,1,window_size);
    window_size = max_region_length/2; %ceil(seq_length * 0.025);
    b=(1/window_size)*ones(1,window_size);

    filtered_signal=filter(b,a,averaged_signal);
    
    plot(filtered_signal,'b','LineWidth',1); 
    set(get(gca, 'XLabel'), 'String', 'Frame number');
    set(get(gca, 'YLabel'), 'String', 'Score');
    hold on;
    
    % define region locations about least similar points
    %local_minima=islocalmin(filtered_signal,'MinProminence',min_prominence);
    local_minima=islocalmin(filtered_signal,'MinSeparation',ceil(max_region_length/2),'SamplePoints',1:length(filtered_signal));
    region_stops=find(local_minima);
    region_stops(end+1) = length(filtered_signal);

    
    % push out last region stop if it's close enough to the end of the
    % sequence
    if region_stops(end) - region_stops(end-1) < max_region_length
        region_stops(end - 1) = region_stops(end);
        region_stops(end) = [];
    end

    msg=sprintf('Found %d regions in %d input frames.\n',length(region_stops), size(descriptors,1));
    disp(msg);
    
    % integrate boundary proximity score
    start_r = 1;
    region_starts = [1];
    for region_id = 1:length(region_stops)
        end_r = region_stops(region_id);
        feas(start_r) = 0;
        feas(start_r+1) = feas(start_r+1) * 0.5;
        feas(end_r) = 0;
        feas(end_r-1) = feas(end_r-1) * 0.5;
        start_r = end_r;
    end
    
    % for each region, set keyframe at most similar frame
    start_r = 1;
    keyframe_indices = int32.empty;
    region_starts = [1];
    for region_id = 1:length(region_stops)
        end_r = region_stops(region_id);
        region = feas(start_r:end_r,:);
        [val,idx] = max(region);
        keyframe_indices(end+1) = start_r + idx - 1; 
        region_starts(end+1) = end_r+1;
        start_r = end_r;
        
        plot(keyframe_indices(end),val,'r.','MarkerSize',40);
        line_h = filtered_signal(start_r);
        line([start_r start_r], [0 line_h],'Color','magenta','LineStyle','--','LineWidth',1);
    end
    
    plot(feas,'r','LineWidth',1); 
    
    % save output to intermediary file
    idcs   = strfind(filepath,'/');
    outdir = filepath(1:idcs(end-1)-1);
    fileID=fopen(outdir +"/input/regions.txt",'w');
    for i = 1:length(keyframe_indices)
        fprintf(fileID,'%d,%d,%d\n',region_starts(i)-1,region_stops(i)-1,keyframe_indices(i)-1);
    end
    fclose(fileID);
    
    leg=legend('Similarity Score','Keyframe Indices','Region Boundaries');
    leg.FontSize = 20;
    ax=gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.75]);
    title('Keyframe Indices and Region Boundaries');
    ax.TitleFontSizeMultiplier = 2;
    saveas(gcf,outdir+"/input/sim_score_filtered.png");
%     
%     hold off
%     plot(filtered_signal,'b','LineWidth',2); 
%     set(get(gca, 'XLabel'), 'String', 'Frame number');
%     set(get(gca, 'YLabel'), 'String', 'Score');
%     hold on
%     feas = load('feas.txt');
%     feas_n = normalize(feas,'range')
%     plot(feas_n,'r','LineWidth',2);
%     leg=legend('Similarity Score','Feasibility Score');
%     leg.FontSize = 20;
%     ax=gca;
%     ax.YAxis.FontSize = 20;
%     ax.XAxis.FontSize = 20;
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.75]);
%     title('Similarity Vs Feasibility')
%     ax.TitleFontSizeMultiplier = 2;
%     saveas(gcf,"sim_feas.png");
    



end
