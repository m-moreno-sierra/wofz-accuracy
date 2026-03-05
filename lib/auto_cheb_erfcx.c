//--- Begin of auto-generated code; do not edit
//
// Generated on 2026-03-04, 11:48:44
// by the piecewise polynomial approximation generator (https://jugit.fz-juelich.de/mlz/ppapp)
// Reference: Wuttke and Kleinsorge,
//            "Code generation for piecewise Chebyshev approximation."
//
// clang-format off
static const double ppapp_a = 0.125; // begin of domain
static const double ppapp_b = 12; // end of domain
static const int ppapp_M = 0; // 2^M subdomains per octave
static const int ppapp_N = 16; // polynomial degree
static const int ppapp_nr = 7; // total number of subdomains
static const int ppapp_j0 = -2; // first octave starts at 2^(j0−1)
static const int ppapp_l0 = 0; // index of a in first octave
static const int ppapp_Nout = 17; // stored coeffs per subdomain
static const int ppapp_nTables = 3;

alignas(64) static const double ppapp_Coeffs0[7 * 8] = {
    0x1.20c8babd2a4cdp-80, -0x1.acd6005a7190dp-76, 0x1.33c0deee255afp-70, -0x1.ae28d7d4df33cp-65, 0x1.21f1403c99ff0p-59, -0x1.78aff2d6a72b2p-54, 0x1.d63a3313cc863p-49, -0x1.1903b154928a7p-43, // subdomain 0:0 (0.125..0.25)
    0x1.9fd48cefce83cp-66, -0x1.3f344721f6a6dp-62, 0x1.d522b6b27bf20p-58, -0x1.555df0c535cbep-53, 0x1.dd5d10c9fd876p-49, -0x1.421e8e7728feep-44, 0x1.a2698a33c9036p-40, -0x1.04aee64535873p-35, // subdomain 1:0 (0.25..0.5)
    0x1.ddf8915290468p-53, -0x1.882a0e5f816c5p-50, 0x1.295f691f09683p-46, -0x1.d89602f5a1955p-43, 0x1.638aa5a57c687p-39, -0x1.02b18506b9467p-35, 0x1.6b9802d95ab58p-32, -0x1.ec0cf48f74a47p-29, // subdomain 2:0 (0.5..1.0)
    0x1.c8f04989af837p-42, -0x1.ab7a43c16d24cp-40, 0x1.4a43dc72e9a64p-37, -0x1.3f0eb2255f02ap-34, 0x1.162a260dd7eacp-31, -0x1.d41d222265876p-29, 0x1.7f947e8f4a7adp-26, -0x1.30c0dae5312c1p-23, // subdomain 3:0 (1.0..2.0)
    0x1.584be0545ee4ep-34, -0x1.9c0d82ba6cb91p-33, 0x1.2e720fef5e461p-31, -0x1.be6ab148d8901p-29, 0x1.05e9dd9a131fap-26, -0x1.1dc0ff719eeb6p-24, 0x1.365bbfc017b0dp-22, -0x1.4b11cdc01882fp-20, // subdomain 4:0 (2.0..4.0)
    0x1.2c552f4beab07p-30, -0x1.0eee34707c913p-29, 0x1.4b777cb7fb114p-29, -0x1.15cf1e63e7c30p-26, 0x1.1507b178c5326p-24, -0x1.c5a4fc13b74cbp-23, 0x1.8169269898684p-21, -0x1.4903e0729cb42p-19, // subdomain 5:0 (4.0..8.0)
    0x1.20b076a3d0205p-29, -0x1.ca747fe7d1fc9p-29, 0x1.c1da407a64159p-30, -0x1.3b2718eb1f07bp-26, 0x1.39ee085e3bdfcp-24, -0x1.bd719f4e3f8a2p-23, 0x1.53ba9de8b503ep-21, -0x1.0a8f8c218b32bp-19, // subdomain 6:0 (8.0..16.0)
};

alignas(64) static const double ppapp_Coeffs1[7 * 8] = {
    0x1.403cb8902f036p-38, -0x1.5a29d8cca9b90p-33, 0x1.60b0a4e35e8f9p-28, -0x1.4ff52d28655f5p-23, 0x1.280377e6506cap-18, -0x1.db731eedadcb8p-14, 0x1.549642cc98b7cp-9, -0x1.a47283b89efb1p-5, // subdomain 0:0 (0.125..0.25)
    0x1.3662c24092b8cp-31, -0x1.5f77477a42cb4p-27, 0x1.78491fa73c187p-23, -0x1.7a16147a55a38p-19, 0x1.611afb945d2dcp-15, -0x1.2e82dbf846fecp-11, 0x1.d28c0e1177cd5p-8, -0x1.3a5c679d7bb59p-4, // subdomain 1:0 (0.25..0.5)
    0x1.3f8189a38a1fap-25, -0x1.8c97dd4edb021p-22, 0x1.d43a7c7a3ecfep-19, -0x1.054d68295b115p-15, 0x1.1192f5bd6881ap-12, -0x1.09e77d40e023ap-9, 0x1.d90093ae10927p-7, -0x1.78cdd551ee51ap-4, // subdomain 2:0 (0.5..1.0)
    0x1.d450afa44eb69p-21, -0x1.5b0ac014164d9p-18, 0x1.ee705e4da7778p-16, -0x1.513ed75fcd78fp-13, 0x1.b65944f3751fap-11, -0x1.0dc51d29420efp-8, 0x1.37ea271bc53d8p-6, -0x1.4f1988444caf6p-4, // subdomain 3:0 (1.0..2.0)
    0x1.5960e2213856ep-18, -0x1.60251dbb01b6dp-16, 0x1.5e738bd670cdcp-14, -0x1.53dec9aae4749p-12, 0x1.409cc2fb85980p-10, -0x1.259061bab4f3ap-8, 0x1.043fe1a9813c1p-6, -0x1.bd6ae4d14af38p-5, // subdomain 4:0 (2.0..4.0)
    0x1.151fe860068aep-17, -0x1.cde5c91445010p-16, 0x1.7d0d9e65876e1p-14, -0x1.370d0385247a5p-12, 0x1.f64cdb5238ec7p-11, -0x1.9106a7d46610dp-9, 0x1.3c7764a6f0203p-7, -0x1.ed7f66d9cb76cp-6, // subdomain 5:0 (4.0..8.0)
    0x1.9dbd2ac7af68cp-18, -0x1.3f8bfcd4864a8p-16, 0x1.ec399c60ec5f9p-15, -0x1.79e915af5e7dbp-13, 0x1.21338dbd03288p-11, -0x1.b92eddd5054cfp-10, 0x1.4f67f46f3b5e3p-8, -0x1.fc477b46bfbd4p-7, // subdomain 6:0 (8.0..16.0)
};

alignas(64) static const double ppapp_Coeffs2[7 * 1] = {
    0x1.a36bbb7f3686dp-1, // subdomain 0:0 (0.125..0.25)
    0x1.5f28ade3ca4acp-1, // subdomain 1:0 (0.25..0.5)
    0x1.038d54ea3d834p-1, // subdomain 2:0 (0.5..1.0)
    0x1.494daffa2ad68p-2, // subdomain 3:0 (1.0..2.0)
    0x1.6e9827d229d58p-3, // subdomain 4:0 (2.0..4.0)
    0x1.7c0348489dbd3p-4, // subdomain 5:0 (4.0..8.0)
    0x1.7fd46c5e09858p-5, // subdomain 6:0 (8.0..16.0)
};
// clang-format on
//--- End of auto-generated code
