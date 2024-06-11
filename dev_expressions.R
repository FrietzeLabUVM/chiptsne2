topfn <- function(df, ex)
{
    eval(substitute(fn(df, ex)))
}

topfn_sub <- function(df, ex)
{
    pr = substitute(ex)
    subset(df, eval(pr))
    # this is how subset handles expressions
    # e <- substitute(ex)
    # r <- eval(e, df, parent.frame())
}

topfn2 <- function(df, ex)
{
    test_expr = substitute(ex)
    ex_str = deparse(test_expr)
    #this worked
    #eval(substitute(fn(df, eval(parse(text = ex_str)))))
    eval(substitute(fn(df, eval(parse(text = ex_str)))))
}


fn <- function(dfr, ex)
{
    eval(substitute(ex), dfr)
}

df <- data.frame(a = 1:5, b = 1:5 )

fn(df, a < 3)
topfn2(df, a < 3)

subset(df, a < 3)
topfn_sub(df, a < 3)



fn(df, a)
topfn(df, a)
topfn2(df, a)

fn(df, 2 * a + b)
topfn(df, 2 * a + b)
topfn2(df, 2 * a + b)




