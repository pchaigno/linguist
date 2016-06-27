module Linguist
  module Strategy
    # Detects language based on filename and/or extension
    class Filename
      def self.call(blob, possible_languages)
        detected_languages = Language.find_by_filename(blob.name.to_s)
        languages = possible_languages & detected_languages
        if languages.empty?
          possible_languages | detected_languages
        else
          languages
        end
      end
    end
  end
end
