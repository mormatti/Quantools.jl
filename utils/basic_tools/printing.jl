"""A function to print colored text in the standard output."""
function colorprint(r, g, b, text)
    print("\e[1m\e[38;2;$r;$g;$b;249m", text)
end

export colorprint