# Make sure the test is always running in the same directory and within the same project.
using Pkg
const MAIN_PROJECT_DIR = joinpath(@__DIR__, "..")
cd(MAIN_PROJECT_DIR)
Pkg.activate(".")

using MarkovJunior; const MJ = MarkovJunior
MJ.markovjunior_asserts_enabled() = true

using Bplus; @using_bplus

const BIG_TEST = @markovjunior 3 'R' begin
    @pragma Hi 1 3 22
    @pragma hello

    @rewrite R => G
    @rewrite RGB => bgw
    @rewrite 3 R => G
    @rewrite (area/2) R => _
    @rewrite (area*2) _ => R
    @rewrite RGB => [2]_[1]

    @rewrite (1.5*area)    R[Gw]B  => b[MR]{wgb}
    @rewrite (length/2.0)  ___     => [1][3][2]  *2
    @rewrite (length*3.0)  [Rw][G] => {wgb}[M]   /1.5

    @rewrite (4.0*length)  R=>G  %0.1
    @rewrite (2:10)        R=>G  %0.2  *3.5
    @rewrite               R=>G  %0.3  /4.1

    @rewrite ((area*4.2):(0.5*length)) RM=>GT  ^[ x ]
    @rewrite (8:(area/4.2))            RM=>GT  ^[ +1 ]
    @rewrite                           RM=>GT                   temperature(0.2)
    @rewrite                           RM=>GT  ^[ x, -2, 4... ] temperature(0.1)

    @rewrite begin
        R => G
        R_[Bb]w => [2]_[bB]{wbR}  %(0.4:0.6)  *0.2   ^[z...]
    end begin
        temperature(40.9)
    end
end
const BIG_TEST_ANSWER = MJ.MarkovAlgorithm(
    MJ.CELL_CODE_BY_CHAR['R'],
    3, 3,

    Pair{Symbol, Vector{Any}}[
        :Hi => Any[ 1, 3, 22 ],
        :hello => Any[ ]
    ],

    [
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            nothing,
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['b']),
                        (MJ.CELL_CODE_BY_CHAR['G'], MJ.CELL_CODE_BY_CHAR['g']),
                        (MJ.CELL_CODE_BY_CHAR['B'], MJ.CELL_CODE_BY_CHAR['w'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            nothing,
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            3,
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.RewriteRuleCell_Wildcard())
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            MJ.ThresholdByArea(0.5f0),
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.RewriteRuleCell_Wildcard(), MJ.CELL_CODE_BY_CHAR['R'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            MJ.ThresholdByArea(2.0f0),
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.RewriteRuleCell_Lookup(2)),
                        (MJ.CELL_CODE_BY_CHAR['G'], MJ.RewriteRuleCell_Wildcard()),
                        (MJ.CELL_CODE_BY_CHAR['B'], MJ.RewriteRuleCell_Lookup(1))
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            nothing,
            ()
        ),

        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['b']),
                        (
                            MJ.RewriteRuleCell_Set('w', 'G'),
                            MJ.RewriteRuleCell_List((
                                MJ.CELL_CODE_BY_CHAR['R'],
                                MJ.CELL_CODE_BY_CHAR['M']
                            ))
                        ),
                        (
                            MJ.CELL_CODE_BY_CHAR['B'],
                            MJ.RewriteRuleCell_Set('w', 'g', 'b')
                        )
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            MJ.ThresholdByArea(1.5f0),
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.RewriteRuleCell_Wildcard(), MJ.RewriteRuleCell_Lookup(1)),
                        (MJ.RewriteRuleCell_Wildcard(), MJ.RewriteRuleCell_Lookup(3)),
                        (MJ.RewriteRuleCell_Wildcard(), MJ.RewriteRuleCell_Lookup(2))
                    ),
                    nothing, 2.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            MJ.ThresholdByLength(0.5f0),
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (
                            MJ.RewriteRuleCell_Set('w', 'R'),
                            MJ.RewriteRuleCell_Set('w', 'g', 'b')
                        ),
                        (
                            MJ.RewriteRuleCell_Set('G'),
                            MJ.RewriteRuleCell_List(tuple(
                                MJ.CELL_CODE_BY_CHAR['M']
                            ))
                        )
                    ),
                    nothing, convert(Float32, 1 / 1.5),
                    MJ.GridDir[ ], 0
                )
            ),
            MJ.ThresholdByLength(3.0f0),
            ()
        ),

        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G'])
                    ),
                    0.1f0, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            MJ.ThresholdByLength(4.0f0),
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G'])
                    ),
                    0.2f0, 3.5f0,
                    MJ.GridDir[ ], 0
                )
            ),
            MJ.ThresholdRange(2, 10),
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G'])
                    ),
                    0.3f0, convert(Float32, 1 / 4.1),
                    MJ.GridDir[ ], 0
                )
            ),
            nothing,
            ()
        ),

        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G']),
                        (MJ.CELL_CODE_BY_CHAR['M'], MJ.CELL_CODE_BY_CHAR['T'])
                    ),
                    nothing, 1.0f0,
                    [ MJ.GridDir(1, -1), MJ.GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdRange(
                MJ.ThresholdByArea(4.2f0),
                MJ.ThresholdByLength(0.5f0)
            ),
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G']),
                        (MJ.CELL_CODE_BY_CHAR['M'], MJ.CELL_CODE_BY_CHAR['T'])
                    ),
                    nothing, 1.0f0,
                    [ MJ.GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdRange(
                8,
                MJ.ThresholdByArea(convert(Float32, 1/4.2))
            ),
            ()
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G']),
                        (MJ.CELL_CODE_BY_CHAR['M'], MJ.CELL_CODE_BY_CHAR['T'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                )
            ),
            nothing,
            tuple(
                MJ.MarkovBiasTemperature(0.2)
            )
        ),
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G']),
                        (MJ.CELL_CODE_BY_CHAR['M'], MJ.CELL_CODE_BY_CHAR['T'])
                    ),
                    nothing, 1.0f0,
                    [ MJ.GridDir(1, -1), MJ.GridDir(1, 1), MJ.GridDir(2, -1) ], 3
                )
            ),
            nothing,
            tuple(
                MJ.MarkovBiasTemperature(0.1)
            )
        ),

        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 0
                ),
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.RewriteRuleCell_Lookup(2)),
                        (MJ.RewriteRuleCell_Wildcard(), MJ.RewriteRuleCell_Wildcard()),
                        (
                            MJ.RewriteRuleCell_Set('b', 'B'),
                            MJ.RewriteRuleCell_List((
                                MJ.CELL_CODE_BY_CHAR['B'],
                                MJ.CELL_CODE_BY_CHAR['b']
                            ))
                        ),
                        (
                            MJ.CELL_CODE_BY_CHAR['w'],
                            MJ.RewriteRuleCell_Set('w', 'b', 'R')
                        )
                    ),
                    (0.4f0, 0.6f0), 0.2f0,
                    MJ.GridDir[ ], 2
                )
            ),
            nothing,
            tuple(
                MJ.MarkovBiasTemperature(40.9f0)
            )
        )
    ]
)

function test_compare(a::MJ.MarkovAlgorithm, b::MJ.MarkovAlgorithm)
    println("Comparing two algorithms...")
    if length(a.sequence) != length(b.sequence)
        println("\tA has ", length(a.sequence), " elements while B has ", length(b.sequence), "!")
    else
        for (oa, ob, i) in zip(a.sequence, b.sequence, 1:length(a.sequence))
            println("\tOp ", i, ", ", typeof(oa), "(", typeof(ob), ")...")
            test_compare(oa, ob)
        end
    end
    return nothing
end
test_compare(a::MJ.AbstractMarkovOp, b::MJ.AbstractMarkovOp) = println("\t\tMismatched/Unsupported types!")
function test_compare(a::MJ.MarkovOpRewrite1D, b::MJ.MarkovOpRewrite1D)
    if length(a.rules) != length(b.rules)
        println("\t\tA has ", length(a.rules), " rules while B has ", length(b.rules), "!")
    else
        for (ra, rb, i) in zip(a.rules, b.rules, 1:length(a.rules))
            println("\t\tRule ", i, ", ", MJ.dsl_string(ra), "    (", MJ.dsl_string(rb), ")")
            if length(ra.cells) != length(rb.cells)
                println("\t\t\tA has ", length(ra.cells), " cells, while B has ", length(rb.cells), "!")
            else
                for (ca, cb, j) in zip(ra.cells, rb.cells, 1:length(ra.cells))
                    println("\t\t\tCell ", j, ", ", ca, "   (", cb, ")")
                    if ca != cb
                        println("\t\t\t\tMismatch!")
                    end
                end
            end
            if ra.weight != rb.weight
                println("\t\t\tWeight mismatch! ", ra.weight, " vs ", rb.weight)
            end
            if ra.mask != rb.mask
                println("\t\t\tMask mismatch! ",
                        typeof(ra.mask), "(", ra.mask, ") vs ",
                        typeof(rb.mask), "(", rb.mask, ")")
            end
            if length(ra.explicit_symmetries) != length(rb.explicit_symmetries)
                println("\t\t\tA has ", length(ra.explicit_symmetries), " explicit symmetries ",
                         "while B has ", length(rb.explicit_symmetries))
            else
                for (sa, sb, j) in zip(ra.explicit_symmetries, rb.explicit_symmetries, 1:length(ra.explicit_symmetries))
                    println("\t\t\tExplicit symmetry ", sa, "  (", sb, ")...")
                    if sa != sb
                        println("\t\t\t\tMismatch!")
                    end
                end
            end
            if ra.unlimited_symmetries_after_axis != rb.unlimited_symmetries_after_axis
                println("\t\t\tInfinite-symmetry mismatch! ",
                            ra.unlimited_symmetries_after_axis,
                            " vs ", rb.unlimited_symmetries_after_axis)
            end

            if ra != rb
                println("\t\t\tFails to match via == operator!")
            end
        end
    end
    if a.threshold != b.threshold
        println("\t\t\tThreshold mismatch! ",
                typeof(a.threshold), "(", a.threshold, ") vs ",
                typeof(b.threshold), "(", b.threshold, ")")
    end
    if length(a.biases) != length(b.biases)
        println("\t\tA has ", length(a.biases), " biases, but B has ", length(b.biases), "!")
    else
        for (ba, bb, i) in zip(a.biases, b.biases, 1:length(a.biases))
            println("\t\tBias ", i, ", ", MJ.dsl_string(ba), "   (", MJ.dsl_string(bb), ")")
            if ba != bb
                println("\t\t\tMismatch!")
            end
        end
        if a.biases != b.biases
            println("\t\tBias tuple Mismatch!")
        end
    end
    if a != b
        println("\t\tFails to match via == operator!")
    end
    return nothing
end

# println("-------------------------------------------------\nEXPECTED:\n",
#         MJ.dsl_string(BIG_TEST_ANSWER), "\n\n\n")
# println("-------------------------------------------------\nACTUAL:\n",
#         MJ.dsl_string(BIG_TEST), "\n\n\n")

@bp_check(BIG_TEST == BIG_TEST_ANSWER, test_compare(BIG_TEST, BIG_TEST_ANSWER))

const BIG_TEST_2 = MJ.parse_markovjunior(MJ.dsl_string(BIG_TEST))
@bp_check(BIG_TEST_2 == BIG_TEST_ANSWER, test_compare(BIG_TEST_2, BIG_TEST_ANSWER))