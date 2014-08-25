%%Show all the characters
function showDictionary(dictionary)

for k=1:size(dictionary);
    subplot(6,8,k),imshow(dictionary{k,2},'InitialMagnification',1000),title(dictionary{k,1});
endfor
endfunction




