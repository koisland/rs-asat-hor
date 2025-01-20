#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Token {
    Number,
    Underscore,
    Hyphen,
    Chimera,
    Other(char),
}

impl From<char> for Token {
    fn from(value: char) -> Self {
        match value {
            '0'..='9' => Token::Number,
            '_' => Token::Underscore,
            '-' => Token::Hyphen,
            '/' => Token::Chimera,
            _ => Token::Other(value),
        }
    }
}
