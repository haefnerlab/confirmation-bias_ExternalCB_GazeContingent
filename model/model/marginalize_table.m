function table = marginalize_table(table,dimensions)
for d=dimensions
    table = sum(table,d);
end
table = squeeze(table);
table = table / sum(table(:));
end