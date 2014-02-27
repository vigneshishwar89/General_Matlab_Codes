function corr_Pitch = corrPitch(pitchContour)
a = load(pitchContour);
time = a(:,1);
pitch = a(:,2);
pitchCent = 1200*log2(pitch./55);
ind0 = find(pitchCent~=-Inf);
pc1 = pitchCent(ind0);
tc1 = time(ind0);
k = 1;
for i = 14:length(pc1)
    med = median(pc1(i-13:i-1));
    d1 = pc1(i)-med;
    if (abs(d1)>=1200)
        pc1(i) = pc1(i)-(1200);
        tc1(i) = tc1(i);
    elseif (abs(d1)>=2400)
        pc1(i) = pc1(i)-(2400);
        tc1(i) = tc1(i);
    elseif (abs(d1)>=3600)
        pc1(i) = pc1(i)-(3600);
        tc1(i) = tc1(i);
    elseif (abs(d1)>=4800)
        pc1(i) = pc1(i)-(4800);
        tc1(i) = tc1(i);
    end
end
%pc1 = medfilt1(pc1,3);
pc1_hz = 55*(2.^(pc1./1200));
plot(time,pitch,'.')
hold on
plot(tc1,pc1_hz,'.r')
corr_Pitch = [tc1 pc1_hz];
dlmwrite('corr_pitch_vignesh_saveri_voice.pitch',corr_Pitch,'delimiter','\t')