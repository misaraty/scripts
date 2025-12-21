function y = xsym(x)
if isnumeric(x)
    x = num2str(x);
end
y = str2sym(x);
end
