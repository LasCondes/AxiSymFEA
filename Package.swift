// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "AxiSymFEA",
    platforms: [
        .macOS(.v13),
        .iOS(.v16)
    ],
    products: [
        .library(
            name: "AxiSymFEA",
            targets: ["AxiSymFEA"]
        ),
    ],
    dependencies: [],
    targets: [
        .target(
            name: "AxiSymFEA",
            dependencies: [],
            linkerSettings: [
                .linkedFramework("Accelerate")
            ]
        ),
        .testTarget(
            name: "AxiSymFEATests",
            dependencies: ["AxiSymFEA"]
        ),
    ]
)
