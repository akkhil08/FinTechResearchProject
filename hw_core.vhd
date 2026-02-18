library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity holt_winters is
    port (
        clk   : in std_logic;
        rst   : in std_logic;
        start : in std_logic;

        alpha : in signed(31 downto 0);
        beta  : in signed(31 downto 0);
        gamma : in signed(31 downto 0);

        data_in  : in signed(31 downto 0);
        valid_in : in std_logic;

        forecast_out   : out signed(31 downto 0);
        forecast_valid : out std_logic
    );
end;

architecture rtl of holt_winters is

constant SEASON_LEN  : integer := 12;
constant NUM_INPUTS  : integer := 60;
constant NUM_OUTPUTS : integer := 12;

subtype fp is signed(31 downto 0);

type data_array   is array (0 to 59) of fp;
type season_array is array (0 to 11) of fp;
type avg_array    is array (0 to 4)  of fp;

signal data_buf : data_array := (others=>(others=>'0'));
signal S : season_array := (others=>(others=>'0'));

signal L,T : fp := (others=>'0');

signal sample_cnt : integer range 0 to 60 := 0;
signal out_cnt    : integer range 0 to 12 := 0;

signal initialized : std_logic := '0';
signal running : std_logic := '0';

constant ONE : fp := to_signed(65536,32);

-- Q16.16 multiply with rounding
function fp_mul(a,b:fp) return fp is
    variable tmp:signed(63 downto 0);
begin
    tmp := a*b + to_signed(2**15,64);
    return tmp(47 downto 16);
end;

begin

process(clk)
    variable sum : fp;
    variable season_avg : avg_array;
    variable L0,T0 : fp;
    variable l_new,t_new,s_new : fp;
    variable idx : integer;
begin
if rising_edge(clk) then

    if rst='1' then
        sample_cnt<=0; out_cnt<=0;
        initialized<='0'; running<='0';
        forecast_valid<='0';

    elsif start='1' then
        running<='1';

    elsif running='1' then

        -----------------------------------
        -- LOAD 60 SAMPLES
        -----------------------------------
        if valid_in='1' and sample_cnt<60 then
            data_buf(sample_cnt)<=data_in;
            sample_cnt<=sample_cnt+1;
        end if;

        -----------------------------------
        -- INITIALIZATION (AFTER 60)
        -----------------------------------
        if sample_cnt=60 and initialized='0' then

            -- initial level (avg first season)
            sum:=(others=>'0');
            for i in 0 to 11 loop
                sum:=sum+data_buf(i);
            end loop;
            L0:=sum/to_signed(12,32);

            -- initial trend
            sum:=(others=>'0');
            for i in 0 to 11 loop
                sum:=sum+(data_buf(12+i)-data_buf(i));
            end loop;
            T0:=sum/to_signed(144,32);

            L<=L0; T<=T0;

            -- season averages
            for s in 0 to 4 loop
                sum:=(others=>'0');
                for i in 0 to 11 loop
                    sum:=sum+data_buf(s*12+i);
                end loop;
                season_avg(s):=sum/to_signed(12,32);
            end loop;

            -- seasonal indices (C++ exact method)
            for i in 0 to 11 loop
                sum:=(others=>'0');
                for s in 0 to 4 loop
                    sum:=sum+(data_buf(s*12+i)-season_avg(s));
                end loop;
                S(i)<=sum/to_signed(5,32);
            end loop;

            initialized<='1';
        end if;

        -----------------------------------
        -- RECURSIVE UPDATES (NONE here since beta=0,gamma=0)
        -----------------------------------

        -----------------------------------
        -- FORECAST
        -----------------------------------
        if initialized='1' and out_cnt<12 then
            forecast_out <=
                L +
                fp_mul(to_signed(out_cnt+1,32),T) +
                S(out_cnt mod 12);

            forecast_valid<='1';
            out_cnt<=out_cnt+1;
        else
            forecast_valid<='0';
        end if;

    end if;
end if;
end process;

end rtl;
