function y = handleSplChars(str,schar)
for k = 1:length(schar)
   ind = findstr(str,schar{k});
   for p = 1:length(ind)
       %ind(p) = ind(p) + p - 1    % handling ever increasing string length
       str = [str(1:ind(p)-1) '_' str((ind(p)+1):end)];
   end
end
y = str;
