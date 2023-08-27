chromosome_coverage = {}

with open("alignment/coverage_per_base.txt", "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        chromosome = fields[0]
        coverage = int(fields[2])

        if chromosome not in chromosome_coverage:
            chromosome_coverage[chromosome] = {"total": coverage, "count": 1}
        else:
            chromosome_coverage[chromosome]["total"] += coverage
            chromosome_coverage[chromosome]["count"] += 1

for chromosome, data in chromosome_coverage.items():
    average_coverage = data["total"] / data["count"]
    print(f"Average coverage for {chromosome}: {average_coverage}")
