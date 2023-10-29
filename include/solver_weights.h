#ifndef SOLVER_WEIGHTS_H
#define SOLVER_WEIGHTS_H

const double rk89_c[] = {
    0.00000000000000000000,
    0.04000000000000000083,
    0.09648736013787360954,
    0.14473104020681040738,
    0.57599999999999995648,
    0.22723265646187659761,
    0.54076734353812339062,
    0.64000000000000001332,
    0.47999999999999998224,
    0.06754000000000000281,
    0.25000000000000000000,
    0.67709201535432428365,
    0.81149999999999999911,
    0.90600000000000002753,
    1.00000000000000000000,
    1.00000000000000000000,
};
const double rk89_b[] = {
    0.01458885278405539637,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00202419788788933252,
    0.21780470845697166848,
    0.12748953408543897692,
    0.22446177454631319192,
    0.17872544912599030997,
    0.07594344758096557846,
    0.12948458791975614446,
    0.02947744761261941737,
    0.00000000000000000000,
};
const double rk89_bh[] = {
    0.02034666655224434684,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00000000000000000000,
    0.00000000000000000000,
    1.06961765098270000784,
    0.07680834711303187456,
    0.11307781868852403995,
    0.25525873579819624570,
    -0.98258980869191636653,
    0.39815458244215140571,
    0.00000000000000000000,
    0.00000000000000000000,
    0.04932600711506839042,
};
const double rk89_a[16][16] = {
    {
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.04000000000000000083,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        -0.01988527319182291017,
        0.11637263332969652319,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.03618276005170260184,
        0.00000000000000000000,
        0.10854828015510781247,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        2.27211426429017748774,
        0.00000000000000000000,
        -8.52688644797639838657,
        6.83077218368622141043,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.05094385535389374386,
        0.00000000000000000000,
        0.00000000000000000000,
        0.17558650498090710990,
        0.00070229612707574678,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.14247836686832848763,
        0.00000000000000000000,
        0.00000000000000000000,
        -0.35417994346686842988,
        0.07595315450295100912,
        0.67651576563371229600,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.07111111111111111105,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.32799092876058982826,
        0.24089796012829906013,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.07124999999999999389,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.32688424515752456667,
        0.11561575484247543777,
        -0.03375000000000000222,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.04822677322465810518,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.03948559980495400246,
        0.10588511619346581416,
        -0.02152006320474309300,
        -0.10453742601833482251,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        -0.02609113435754923521,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.03333333333333333287,
        -0.16525040066381049830,
        0.03434664118368616764,
        0.15957582832152089614,
        0.21408573218281934381,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        -0.03628423396255658906,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        -1.09616759742720870641,
        0.18260355043213311044,
        0.07082254444170683894,
        -0.02313647018482431136,
        0.27112047263209326786,
        1.30813374942298077386,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        -0.50746350564169750985,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        -6.63134219865723739673,
        -0.25274801009088010417,
        -0.49526123800360954963,
        0.29325255452538867562,
        1.44010869376828098964,
        6.23793449864705618069,
        0.72701920545269871354,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.61301182569559320434,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        9.08880389164046320616,
        -0.40737881562934485924,
        1.79073338949037474954,
        0.71492716676175505075,
        -1.43858085784172295973,
        -8.26332931206474086139,
        -1.53757057080886516687,
        0.34538328275648716437,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        -1.21169791034387386297,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        -19.05581871559595441568,
        1.26306067538987520926,
        -6.91391696917845788306,
        -0.67646226650949803361,
        3.36786044502660786293,
        18.00675164312590936788,
        6.83882892679427989435,
        -1.03151645192195040579,
        0.41291062321306226668,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        2.15738900749405360102,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        23.80712219809580432184,
        0.88627792492165557992,
        13.13913039759876433266,
        -2.60441570928771470861,
        -5.19385994978387266485,
        -20.41234071154150697680,
        -12.30085625250572256562,
        1.52155309500853941351,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
};

const double rk56_c[] = {
    0.00000000000000000000,
    0.17999999999999999334,
    0.16666666666666665741,
    0.25000000000000000000,
    0.53000000000000002665,
    0.59999999999999997780,
    0.80000000000000004441,
    1.00000000000000000000,
    1.00000000000000000000,
};
const double rk56_b[] = {
    0.07638888888888889506,
    0.00000000000000000000,
    0.00000000000000000000,
    0.36940836940836940805,
    0.00000000000000000000,
    0.24801587301587302292,
    0.23674242424242425420,
    0.06944444444444444753,
    0.00000000000000000000,
};
const double rk56_bh[] = {
    0.05870020964360587318,
    0.00000000000000000000,
    0.00000000000000000000,
    0.48072562358276643701,
    -0.85341242076919088255,
    1.20464852607709760335,
    0.00000000000000000000,
    -0.05924237307216030646,
    0.16858043453788135180,
};
const double rk56_a[9][9] = {
    {
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.17999999999999999334,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.08950617283950616787,
        0.07716049382716048954,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.06250000000000000000,
        0.00000000000000000000,
        0.18750000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.31651600000000001955,
        0.00000000000000000000,
        -1.04494799999999998796,
        1.25843199999999999505,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.27232612736485628524,
        0.00000000000000000000,
        -0.82513360323886641989,
        1.04809176788124158719,
        0.10471570799276856689,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        -0.16699418599716514544,
        0.00000000000000000000,
        0.63170850202429151832,
        0.17461044552773877236,
        -1.06653564590860661099,
        1.22721088435374148240,
        0.00000000000000000000,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.36423751686909583425,
        0.00000000000000000000,
        -0.20404858299595141080,
        -0.34883737816068643989,
        3.26193230328568661847,
        -2.75510204081632670281,
        0.68181818181818176772,
        0.00000000000000000000,
        0.00000000000000000000,
    },
    {
        0.07638888888888889506,
        0.00000000000000000000,
        0.00000000000000000000,
        0.36940836940836940805,
        0.00000000000000000000,
        0.24801587301587302292,
        0.23674242424242425420,
        0.06944444444444444753,
        0.00000000000000000000,
    },
};

#endif