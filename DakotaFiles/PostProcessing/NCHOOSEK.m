function evals = NCHOOSEK(a,b)
% Vectorized form of MatLab's "nchoosek(a,b)"
evals = factorial(a)./(factorial(a-b).*factorial(b));
end