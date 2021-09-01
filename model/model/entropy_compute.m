function H  = entropy_compute(table,dimensions)
table = table/sum(table(:));
d = 1:ndims(table);
d = setdiff(d,dimensions);
table = marginalize_table(table,d);
H = -dot(table(:), log(table(:)));
end