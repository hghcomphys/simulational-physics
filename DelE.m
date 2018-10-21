function f = DelE(i,j)

global S L;

if i~=1
    im = i-1;
else
    im = L;
end

if j~=1
    jm = j-1;
else
    jm = L;
end

f   = 2 * S(i,j) * ( S(im,j) + S(mod(i,L)+1,j) + S(i,jm) + S(i,mod(j,L)+1) );
