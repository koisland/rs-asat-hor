#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Token {
    SF,
    Chrom,
    Monomer,
    Divergent,
    Live,
    MType,
    Number,
    Chimera,
    Hyphen,
    Value(char),
}

impl Token {
    pub fn char(&self) -> char {
        match self {
            Token::SF => 'S',
            Token::Chrom => 'C',
            Token::Number => 'n',
            Token::Divergent => 'd',
            Token::Live => 'L',
            Token::MType => 'H',
            Token::Monomer => '.',
            Token::Chimera => '/',
            Token::Hyphen => '-',
            Token::Value(value) => *value,
        }
    }
}

impl From<char> for Token {
    fn from(value: char) -> Self {
        match value {
            'S' => Token::SF,
            'C' => Token::Chrom,
            'H' => Token::MType,
            '0'..='9' => Token::Number,
            'd' => Token::Divergent,
            'L' => Token::Live,
            '.' => Token::Monomer,
            '/' => Token::Chimera,
            '-' => Token::Hyphen,
            _ => Token::Value(value),
        }
    }
}
