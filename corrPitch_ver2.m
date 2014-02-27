function corr_Pitch = corrPitch_ver2(pitchContour)
a = load(pitchContour);
time = a(:,1);
pitch = a(:,2);
ind0 = find(pitch~=0);
% pitchCent = 1200*log2(pitch./55);
% ind0 = find(pitchCent~=-Inf);
pc1 = pitch(ind0);
tc1 = time(ind0);
pc2 = pc1;
tc2 = tc1;
% tc1 = time(ind0);
% tc2 = tc1;
k = 1;
for l = 1:1
    for i = 1:length(pc1)-1
        rat1 = pc1(i+1)/pc2(i);
        if (abs(rat1-2)<0.1)
            pc2(i+1) = pc2(i+1)/2;
        elseif (abs(rat1-3)<0.1)
            pc2(i+1) = pc2(i+1)/3;
        elseif (abs(rat1-4)<0.1)
            pc2(i+1) = pc2(i+1)/4;
        elseif (abs(rat1-0.5)<0.05)
            pc2(i+1) = pc2(i+1)*2;
        end
    end
end
%pc1 = medfilt1(pc1,3);
% pc2_hz = 55*(2.^(pc2./1200));
plot(time,pitch,'.')
hold on
plot(tc2,pc2,'r')
corr_Pitch = [tc2 pc2];
dlmwrite('corr_pitch_vignesh_saveri_voice.pitch',corr_Pitch,'delimiter','\t')