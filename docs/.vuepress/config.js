module.exports = {
  title: "Variant calling for NGS data",
  description:
    "A tutorial on variant calling for NGS data using FreeBayes and Delly",
  base: "/courses/vc/",
  themeConfig: {
    repo: "tobiasrausch/variant-calling",
    nav: [
      { text: "Home", link: "/" },
      { text: "SNV Calling", link: "/snv/" },
      { text: "SV Calling", link: "/sv/" },
      { text: "Appendix", link: "/appendix/" }
    ],
    sidebar: ["/snv/", "/sv/", "/appendix/"]
  },
  plugins: {
    "@vuepress/back-to-top": true
  }
};
