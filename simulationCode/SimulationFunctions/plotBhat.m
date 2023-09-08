%% Heatmap of Bhat

function out = plotBhat(type, ss, i, p, freq)

file  = strcat(type(1), '_Bhat_ss', num2str(i), '.txt');
Bhats = dlmread(file);

if strcmp(freq,'low') 
    jj_seq = 25;
else
    jj_seq = [25 75];
end

for jj = jj_seq
    Bhat = Bhats(((jj-1)*p+1):jj*p,:);
    MyHeatmapRed(Bhat);
    figuretitle = strcat('Bhat from', {' '}, type, ': s =', {' '}, num2str(ss(i)));
    title(figuretitle);  
    figurename  = strcat(type(1), '_Bhat_ss', num2str(i), '_', num2str(jj), '.pdf');
    print(figurename,'-dpdf', '-fillpage');
end

end