# Contribute to seapodym-codebase

## Clang-Format

Compare styles [here](https://zed0.co.uk/clang-format-configurator/).

We are using `Google` style with :

- `IndentWidth: 4`
- `AlignAfterOpenBracket: AlwaysBreak`

Apply on all files :

```bash
find . | grep -i ".*\.\(cpp\|hpp\|h\|hpp\)$" | xargs clang-format --style=file -i
````
