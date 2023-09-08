%% Heatmap of predErr

function out = plotCVmesh(type, ss, i, freq)

file  = strcat(type(1), '_predErr_ss', num2str(i), '.txt');
cvms  = dlmread(file);

if strcmp(freq,'low') 
    jj_seq = 25;
else
    jj_seq = [25 75];
end

for jj = jj_seq
    cvm = cvms(((jj-1)*15+1):jj*15,:);
    MyHeatmapRed(cvm);
    hold on;
    [MM, II] = min(cvm);
    Xindex   = find(MM == min(MM));
    Yindex   = II(Xindex);
    plot([Xindex-0.5 Xindex+0.5],[Yindex-0.5 Yindex-0.5], 'Color', [0, 51, 204]/255, 'LineWidth', 2)
    plot([Xindex-0.5 Xindex+0.5],[Yindex+0.5 Yindex+0.5], 'Color', [0, 51, 204]/255, 'LineWidth', 2)
    plot([Xindex-0.5 Xindex-0.5],[Yindex-0.5 Yindex+0.5], 'Color', [0, 51, 204]/255, 'LineWidth', 2)
    plot([Xindex+0.5 Xindex+0.5],[Yindex-0.5 Yindex+0.5], 'Color', [0, 51, 204]/255, 'LineWidth', 2)
    figuretitle = strcat('CV mesh from', {' '}, type, ': s =', {' '}, num2str(ss(i)));
    title(figuretitle);  
    figurename  = strcat(type(1), '_predErr_ss', num2str(i), '_', num2str(jj), '.pdf');
    print(figurename,'-dpdf', '-fillpage');
end

end