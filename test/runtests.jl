# Make sure the test is always running in the same directory and within the same project.
using Pkg
const MAIN_PROJECT_DIR = joinpath(@__DIR__, "..")
cd(MAIN_PROJECT_DIR)
Pkg.activate(".")

using MarkovJunior; const MJ = MarkovJunior
using Bplus; @using_bplus

const BIG_TEST = @markovjunior 3 'R' begin
    @pragma Hi 1 3 22
    @pragma hello

    @rewrite R => G
    @rewrite RGB => bgw
    @rewrite 3 R => G
    @rewrite (area/2) R => _
    @rewrite (area*2) _ => R

    @rewrite (1.5*area)    R[Gw]B  => b[MR]{wgb}
    @rewrite (length/2.0)  ___     => [1][3][2]  *2
    @rewrite (length*3.0)  [Rw][G] => {wgb}[M]   /1.5

    @rewrite (4.0*length)  R=>G  %0.1
    @rewrite (2:10)        R=>G  %0.2  *3.5
    @rewrite               R=>G  %0.3  /4.1

    @rewrite ((area*4.2):(0.5*length)) RM=>GT  S[ x ]
    @rewrite (8:(area/4.2))            RM=>GT  S[ +1 ]
    @rewrite                           RM=>GT  S[ x, -2, 4... ]

    #TODO: Rewrite groups
    #TODO: Biases
end

@bp_check(BIG_TEST == MJ.MarkovAlgorithm(
    MJ.CELL_CODE_BY_CHAR['R'],
    3, 3,

    [
        :Hi => [ 1, 3, 22 ],
        :hello => [ ]
    ],

    [
        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['G'])
                    ),
                    nothing, 1.0f0,
                    [ ], 0
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
                    [ ], 0
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
                    [ ], 0
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
                    [ ], 0
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
                    [ ], 0
                )
            ),
            MJ.ThresholdByArea(2.0f0),
            ()
        ),

        MJ.MarkovOpRewrite1D(
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['b']),
                        (
                            MJ.RewriteRuleCell_List(
                                MJ.CELL_CODE_BY_CHAR['G'],
                                MJ.CELL_CODE_BY_CHAR['w']
                            ),
                            MJ.RewriteRuleCell_List(
                                MJ.CELL_CODE_BY_CHAR['M'],
                                MJ.CELL_CODE_BY_CHAR['R']
                            )
                        ),
                        (
                            MJ.CELL_CODE_BY_CHAR['B'],
                            MJ.RewriteRuleCell_Set('w', 'g', 'b')
                        )
                    ),
                    nothing, 1.0f0,
                    [ ], 0
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
                    [ ], 0
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
                            MJ.RewriteRuleCell_List(
                                MJ.CELL_CODE_BY_CHAR['R'], MJ.CELL_CODE_BY_CHAR['w']
                            ),
                            MJ.RewriteRuleCell_Set('w', 'g', 'b')
                        ),
                        (
                            MJ.RewriteRuleCell_List(
                                MJ.CELL_CODE_BY_CHAR['G']
                            ),
                            MJ.RewriteRuleCell_List(
                                MJ.CELL_CODE_BY_CHAR['M']
                            )
                        )
                    ),
                    nothing, 1.0f0/1.5f0,
                    [ ], 0
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
                    0.1, 1.0f0,
                    [ ], 0
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
                    0.2, 3.5f0,
                    [ ], 0
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
                    0.3, 1.0f0 / 4.1f0,
                    [ ], 0
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
                    [ GridDir(1, -1), GridDir(1, 1) ], nothing
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
                    [ GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdRange(
                8,
                MJ.ThresholdByArea(1.0f0 / 4.2f0)
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
                    [ GridDir(1, -1), GridDir(1, 1), GridDir(2, -1) ], 3
                )
            ),
            nothing,
            ()
        )
    ]
))

@bp_check(MJ.parse_markovjunior(MJ.dsl_string(BIG_TEST)) == BIG_TEST)